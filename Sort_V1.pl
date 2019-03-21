#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(basename);
use List::MoreUtils qw(uniq none);
use POSIX qw(ceil);

=head1 Usage
 options:
  -k                length of kmer, defaults is 17
  -jellyfish   *    files of kmers from jellyfish
  -filelist    *    list of sequence data file(s) in the absolute path
  -group            number of groups you want to sort, defaults is 5
  -length      *    length of reads
  -iteration        number of iteration times, defaults is 1
  -help             output help information to screen
 Example:
  perl Sort_V1.pl -k 17 -jellyfish fit_kmer.txt -filelist list_of_fastq_files -group 100 -length 125 -iteration 1
=cut

my($K,$Jellyfish,$File,$Group,$Number,$Length,$Iteration,$Help);
GetOptions("k:i"         => \$K,             # kmer����
		   "jellyfish:s" => \$Jellyfish,     # ��jellyfish�õ���kmers
		   "filelist:s"  => \$File,          # fastq�б��ļ�
		   "group:i"     => \$Group,         # ������Ŀ
		   "length:i"    => \$Length,        # �趨����fq�ļ���reads�ĳ��ȣ����ڼ���$Number
		   "iteration:i" => \$Iteration,     # �趨��������
		   "help|?"      => \$Help,)
or die("Error in command line arguments!\n");
die `pod2text $0` if($Help || !$File || !$Jellyfish || !$Length);  # δ��������fastq�ļ��б��򱨴�

$K||=17;                                               # Ĭ�ϵ�kmer����Ϊ17
$Group||=5;                                            # Ĭ�ϵķ�����ĿΪ5
$Iteration||=1;                                        # Ĭ�ϵ���һ�Σ���������һ��
$Number=ceil(($Length-$K+1)*0.2);                      # Ĭ�ϵ�read��kmer���ص�������ֵΪ22((125-17+1))*0.2=21.8,����ȡ��Ϊ22
$|=1;                                                  # ��print��write���֮��ǿ����ջ�����

# ===========
# ��jellyfish�õ���uniqe kmers���뵽��ϣ����Ϊkey��valueΪ����������
# ===========
if($Jellyfish=~/\.gz$/){
	open IN,"zcat $Jellyfish |" or die "$!";
}else{
	open IN,"<$Jellyfish" or die "$!";
}

my $jellyfish;
while(<IN>){
	chomp;
	next if /^\s+$/;
	$jellyfish->{$_}=[];
}
print "Kmers from jellyfish have been read in!\n";
close IN;

# ===========
# ��ʼ��
# ===========
open FL, "<$File" or die "$!";                         # ��ȡfastq�ļ���
LINE:while(my $fq = <FL>){
	chomp($fq);
	next if $fq=~/^\s+$/;
	next if $fq=~/(.*?)2\.(fastq|fq)\.gz$/;
	next if $fq=~/(.*?)2\.(fastq|fq)$/;

	if($fq=~/(.*?)1\.(fastq|fq)\.gz$/){                # �򿪳ɶԵ�fastq�ļ�
		my $fq1 = $fq;
		my $fq2 = $1."2.".$2.".gz";
		open FQ1, "zcat $fq1 |" or die "$!";
		open FQ2, "zcat $fq2 |" or die "$!";
	}elsif($fq=~/(.*?)1\.(fastq|fq)$/){
		my $fq1 = $fq;
		my $fq2 = $1."2.".$2;
		open FQ1, $fq1 or die "$!";
		open FQ2, $fq2 or die "$!";
	}else{
		warn "$fq isn't end with .(fastq|fq).gz or .(fastq|fq), maybe it's wrong file!\n";
	}
	
	my $set=1;
	my $count=0;
	my @kmers_fit;
	while ( my $header1 = <FQ1> ) {                         # ��fq1����ȡ4��
		my $seq1 = <FQ1>;
		my $qual_header1 = <FQ1>;
		my $qual1 = <FQ1>;
		die "ERROR: expected \'\@\' but saw $header1" unless substr( $header1, 0, 1 ) eq '@';
		
		my $header2 = <FQ2>;                                # ��fq2����ȡ4��
		my $seq2 = <FQ2>;
		my $qual_header2 = <FQ2>;
		my $qual2 = <FQ2>;
		die "ERROR: expected \'\@\' but saw $header2" unless substr( $header2, 0, 1 ) eq '@';
		
		chomp(my $seq1_no_linefeed = $seq1); 			    # ȥ�����з�
		for(my $i = 0;$i <= length($seq1_no_linefeed)-$K;$i++){
			my $kmer = substr($seq1_no_linefeed,$i,$K);
			
			my $revcom = reverse $kmer;
			$revcom =~ tr/ACGT/TGCA/;
			
			if(exists $jellyfish->{$kmer}){                 # Kmer������������fit_kmer�о���1
				push @kmers_fit,$kmer;
				$count++;
			}elsif(exists $jellyfish->{$revcom}){
				push @kmers_fit,$revcom;
				$count++;
			}
		}
		
		chomp(my $seq2_no_linefeed = $seq2); 			    # ȥ�����з�
		for(my $i = 0;$i <= length($seq2_no_linefeed)-$K;$i++){
			my $kmer = substr($seq2_no_linefeed,$i,$K);
			
			my $revcom = reverse $kmer;
			$revcom =~ tr/ACGT/TGCA/;
			
			if(exists $jellyfish->{$kmer}){
				push @kmers_fit,$kmer;                      # ��һ��reads������uniq kmers��������ͬһ����            
				$count++;
			}elsif(exists $jellyfish->{$revcom}){
				push @kmers_fit,$revcom;
				$count++;
			}
		}

		if($count > ($Length-$K+1)*2*0.8){      # Ѱ��������(125-17+1)*2*0.8=174.4��kmerΪuniq kmer��reads���г�ʼ��,2����һ��reads
			for my $kmer (@kmers_fit){          # ���ö�reads���е�uniq kmers��ֵ�����г�ʼ��
				push @{$jellyfish->{$kmer}},$set;
				@{$jellyfish->{$kmer}} = uniq @{$jellyfish->{$kmer}};
			}
			$set++;			
		}
		$count=0;
		@kmers_fit=();
		last LINE if $set > $Group;                         # �ҵ�Group�����˳�ѭ��
	}
}
print "Initialization is finished!\n";
close FL;
close FQ1;
close FQ2;

open OUT,">","Pool_initialization.txt";
for(keys %$jellyfish){
	if(@{$jellyfish->{$_}}){
		print OUT "$_\t", "@{$jellyfish->{$_}}\n";
	}
}
close OUT;

# ===================
# ��ʼ����Ϻ����¶�ȡ�ļ�
# ===================
my(@kmers,%dataset1,%dataset2,%seq_fq,$fq1,$fq2);
open FL, "<$File" or die "$!";                          # ��ȡfastq�ļ��б�
while(my $fq = <FL>){
	chomp($fq);
	next if $fq=~/^\s+$/;
	next if $fq=~/(.*?)2\.(fastq|fq)\.gz$/;
	next if $fq=~/(.*?)2\.(fastq|fq)$/;
	
	if($fq=~/(.*?)1\.(fastq|fq)\.gz$/){                 # �򿪳ɶԵ�fastq�ļ�
		$fq1 = $fq;
		$fq2 = $1."2.".$2.".gz";
		open FQ1, "zcat $fq1 |" or die "$!";
		open FQ2, "zcat $fq2 |" or die "$!";
	}elsif($fq=~/(.*?)1\.(fastq|fq)$/){
		$fq1 = $fq;
		$fq2 = $1."2.".$2;
		open FQ1, $fq1 or die "$!";
		open FQ2, $fq2 or die "$!";
	}else{
		warn "$fq isn't end with .(fastq|fq).gz or .(fastq|fq), maybe it's wrong file!\n";
	}	
	
	open REMAIN,">>Remain1.fq" or die "$!";
	
	while ( my $header1 = <FQ1> ) {                                      # ��fq1����ȡ4��
		my $seq1 = <FQ1>;
		my $qual_header1 = <FQ1>;
		my $qual1 = <FQ1>;
		#die "ERROR: expected \'\@\' but saw $header1" unless substr( $header1, 0, 1 ) eq '@';
		
		my $header2 = <FQ2>;                                             # ��fq2����ȡ4��
		my $seq2 = <FQ2>;
		my $qual_header2 = <FQ2>;
		my $qual2 = <FQ2>;
		#die "ERROR: expected \'\@\' but saw $header2" unless substr( $header2, 0, 1 ) eq '@';

		chomp(my $seq1_no_linefeed=$seq1);                               # ȥ��seq����Ļ��з�
		$seq_fq{$seq1_no_linefeed}=$header1.$seq1.$qual_header1.$qual1;  # ��fastq�ļ�ת��Ϊ����Ϊ��������������fastq�ļ�Ϊֵ��hash��
		
		chomp(my $seq2_no_linefeed=$seq2);
		$seq_fq{$seq2_no_linefeed}=$header2.$seq2.$qual_header2.$qual2;

		for(my $i=0; $i <= length($seq1_no_linefeed)-$K; $i++){
			my $kmer = substr($seq1_no_linefeed,$i,$K);
			
			my $revcom = reverse $kmer;
			$revcom =~ tr/ACGT/TGCA/;
			
			if(exists $jellyfish->{$kmer}){
				for my $num_set (@{$jellyfish->{$kmer}}){
					$dataset1{$num_set}++;                       # ���read��ĳ��kmerλ��kmer���У������ݼ��ۼơ�
				}
				push @kmers,$kmer;
			}elsif(exists $jellyfish->{$revcom}){
				for my $num_set (@{$jellyfish->{$revcom}}){
					$dataset1{$num_set}++;
				}
				push @kmers,$revcom;		
			}
		}

		for(my $i=0; $i <= length($seq2_no_linefeed)-$K; $i++){
			my $kmer = substr($seq2_no_linefeed,$i,$K);
			
			my $revcom = reverse $kmer;
			$revcom =~ tr/ACGT/TGCA/;
			
			if(exists $jellyfish->{$kmer}){
				for my $num_set (@{$jellyfish->{$kmer}}){
					$dataset2{$num_set}++;                       # ���read��ĳ��kmerλ��kmer���У������ݼ��ۼơ�
				}
				push @kmers,$kmer;
			}elsif(exists $jellyfish->{$revcom}){
				for my $num_set (@{$jellyfish->{$revcom}}){
					$dataset2{$num_set}++;
				}
				push @kmers,$revcom;
			}
		}

		my @fit1=grep{$dataset1{$_}>$Number}(keys %dataset1);
		my @fit2=grep{$dataset2{$_}>$Number}(keys %dataset2);

		my @fit=uniq @fit1,@fit2;                               # һ��reads�У�ֻҪ��һ��������������

		if($#fit==-1){                                          # ���������������remain���ݼ���
			print REMAIN $seq_fq{$seq1_no_linefeed};
			print REMAIN $seq_fq{$seq2_no_linefeed};
		}
		elsif($#fit==0){                                        # ���ֻ��һ��Ԫ�����㣬����read������kmers�����kmer��
			my $file="Group".$fit[0].".fq"; 
			open OUT, ">>$file" or die "$!";
			print OUT $seq_fq{$seq1_no_linefeed};
			print OUT $seq_fq{$seq2_no_linefeed};
			for my $kmer (@kmers){
				push @{$jellyfish->{$kmer}},$fit[0];
				@{$jellyfish->{$kmer}} = uniq @{$jellyfish->{$kmer}};
			}
		}
		else{                                                   # ����ж��Ԫ������
			my $first=shift(@fit);                              # �޸�kmer�������ݼ�ָ���һ���������������ݼ�
			my $file="Group".$first.".fq";                      # ȡ��һ����������ݼ��ı��
			
			open OUT, ">>$file" or die "$!";
			print OUT $seq_fq{$seq1_no_linefeed};
			print OUT $seq_fq{$seq2_no_linefeed};
			
			unless(system "for i in @fit; do cat Group\$i.fq >> $file; done"){
				for(@fit){
					unlink "Group$_.fq" or warn "$!";           # ���гɹ���ɾ��ԭ�ļ� 
					print "Group$_.fq is merged with $file and Group$_.fq has been deleted successfully!\n";
				}
			}

			for my $kmer (@kmers){
				push @{$jellyfish->{$kmer}},$first;
				@{$jellyfish->{$kmer}} = uniq @{$jellyfish->{$kmer}};
			}
			
			for my $num_set (@fit){                             # �ж���������������ֵ,��ȥ��
				for my $kmer (keys %$jellyfish){
					for(@{$jellyfish->{$kmer}}){
						$_ = $first if $_ == $num_set;
					}
					@{$jellyfish->{$kmer}} = uniq @{$jellyfish->{$kmer}};
				}
			}
		}
		(%seq_fq,%dataset1,%dataset2,@fit,@fit1,@fit2,@kmers)=();                     # ����Խ�ʡ�ڴ�
	}
	print "Round 1 iteration of ", basename($fq1)," is finished!\n";
	print "Round 1 iteration of ", basename($fq2)," is finished!\n";
}
close FQ1;
close FQ2;
close OUT;
close REMAIN;
close FL;

open OUT,">","Pool_round1.txt";
for(keys %$jellyfish){
	if(@{$jellyfish->{$_}}){
		print OUT "$_\t", "@{$jellyfish->{$_}}\n";
	}
}
close OUT;

#==================
# ��������Remain.fq��Ѱ��
#==================
my %dataset;
for(my $j=2;$j<=$Iteration;$j++){
	my $p=$j-1;
	open FQ,"<Remain$p.fq" or die "$!";                         # ��Remain1.fq�е�reads������Ѱ
	open REMAIN,">>Remain$j.fq" or die "$!";                    # ��ʣ�µ�reads�����Remain2.fq��
	while ( my $header = <FQ> ) {
		my $seq = <FQ>;
		my $qual_header = <FQ>;
		my $qual = <FQ>;
		#die "ERROR: expected \'\@\' but saw $header" unless substr( $header, 0, 1 ) eq '@';
	
		chomp(my $seq_no_linefeed=$seq);                               # �˴�ȥ��seq����Ļ��з�
		$seq_fq{$seq_no_linefeed}=$header.$seq.$qual_header.$qual;     # ��fastq�ļ�ת��Ϊ����Ϊ��������������fastq�ļ�Ϊֵ��hash��
	
		for(my $i=0;$i<=length($seq_no_linefeed)-$K;$i++){
			my $kmer=substr($seq_no_linefeed,$i,$K);
			
			my $revcom = reverse $kmer;
			$revcom =~ tr/ACGT/TGCA/;
			
			if(exists $jellyfish->{$kmer}){
				for my $num_set (@{$jellyfish->{$kmer}}){
					$dataset{$num_set}++;                       # ���read��ĳ��kmerλ��kmer���У������ݼ��ۼơ�
				}
				push @kmers,$kmer;
			}elsif(exists $jellyfish->{$revcom}){
				for my $num_set (@{$jellyfish->{$revcom}}){
					$dataset{$num_set}++;
				}
				push @kmers,$revcom;
			}
		}

		my @fit=grep{$dataset{$_}>$Number}(keys %dataset);      # ���������ֵ���򽫸����ݼ���ż���@fit��
		
		if($#fit==-1){                                          # ���������������remain���ݼ���
			print REMAIN $seq_fq{$seq_no_linefeed};
		}
		elsif($#fit==0){                                        # ���ֻ��һ��Ԫ������
			my $file="Group".$fit[0].".fq";                     # ȡ�����������ݼ��ı��
			open OUT, ">>$file" or die "$!";
			print OUT $seq_fq{$seq_no_linefeed};
			for my $kmer (@kmers){
				push @{$jellyfish->{$kmer}},$fit[0];
				@{$jellyfish->{$kmer}} =  uniq @{$jellyfish->{$kmer}};
			}
		}
		else{                                                   # ����ж��Ԫ������
			my $first=shift(@fit);                              # �޸�kmer�������ݼ�ָ���һ���������������ݼ�
			my $file="Group".$first.".fq";                      # ȡ��һ����������ݼ��ı��
			
			open OUT, ">>$file" or die "$!";
			print OUT $seq_fq{$seq_no_linefeed};
			
			unless(system "for i in @fit; do cat Group\$i.fq >> $file; done"){
				for(@fit){
					unlink "Group$_.fq" or warn "$!";           # ���гɹ���ɾ��ԭ�ļ� 
					print "Group$_.fq is merged with $file and Group$_.fq has been deleted successfully!\n";
				}
			}

			for my $kmer (@kmers){
				push @{$jellyfish->{$kmer}},$first;
				@{$jellyfish->{$kmer}} = uniq @{$jellyfish->{$kmer}};
			}
			
			for my $num_set (@fit){                             # �ж���������������ֵ,��ȥ��
				for my $kmer (keys %$jellyfish){
					for(@{$jellyfish->{$kmer}}){
						$_ = $first if $_ == $num_set;
					}
					@{$jellyfish->{$kmer}} = uniq @{$jellyfish->{$kmer}};
				}
			}
		}
		(%seq_fq,%dataset,@fit,@kmers)=();                     # ����Խ�ʡ�ڴ�		
	}
	close FQ;
	close OUT;
	close REMAIN;
	
	unlink "Remain$p.fq" or warn "$!";                         # ɾ���ϴε�����Remain.fq
	print "Round $j iteration is finished!\n";
}

open OUT,">","Pool_finish.txt";
for(keys %$jellyfish){
	if(@{$jellyfish->{$_}}){
		print OUT "$_\t", "@{$jellyfish->{$_}}\n";
	}
}
close OUT;
__END__
