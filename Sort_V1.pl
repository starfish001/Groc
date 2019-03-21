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
GetOptions("k:i"         => \$K,             # kmer长度
		   "jellyfish:s" => \$Jellyfish,     # 从jellyfish得到的kmers
		   "filelist:s"  => \$File,          # fastq列表文件
		   "group:i"     => \$Group,         # 分类数目
		   "length:i"    => \$Length,        # 设定输入fq文件中reads的长度，用于计算$Number
		   "iteration:i" => \$Iteration,     # 设定迭代次数
		   "help|?"      => \$Help,)
or die("Error in command line arguments!\n");
die `pod2text $0` if($Help || !$File || !$Jellyfish || !$Length);  # 未定义输入fastq文件列表则报错

$K||=17;                                               # 默认的kmer长度为17
$Group||=5;                                            # 默认的分类数目为5
$Iteration||=1;                                        # 默认迭代一次，即仅运行一次
$Number=ceil(($Length-$K+1)*0.2);                      # 默认的read与kmer池重叠条数阈值为22((125-17+1))*0.2=21.8,向上取整为22
$|=1;                                                  # 在print和write语句之后强制清空缓存区

# ===========
# 将jellyfish得到的uniqe kmers读入到哈希中作为key，value为空匿名数组
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
# 初始化
# ===========
open FL, "<$File" or die "$!";                         # 读取fastq文件列
LINE:while(my $fq = <FL>){
	chomp($fq);
	next if $fq=~/^\s+$/;
	next if $fq=~/(.*?)2\.(fastq|fq)\.gz$/;
	next if $fq=~/(.*?)2\.(fastq|fq)$/;

	if($fq=~/(.*?)1\.(fastq|fq)\.gz$/){                # 打开成对的fastq文件
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
	while ( my $header1 = <FQ1> ) {                         # 打开fq1，读取4行
		my $seq1 = <FQ1>;
		my $qual_header1 = <FQ1>;
		my $qual1 = <FQ1>;
		die "ERROR: expected \'\@\' but saw $header1" unless substr( $header1, 0, 1 ) eq '@';
		
		my $header2 = <FQ2>;                                # 打开fq2，读取4行
		my $seq2 = <FQ2>;
		my $qual_header2 = <FQ2>;
		my $qual2 = <FQ2>;
		die "ERROR: expected \'\@\' but saw $header2" unless substr( $header2, 0, 1 ) eq '@';
		
		chomp(my $seq1_no_linefeed = $seq1); 			    # 去掉换行符
		for(my $i = 0;$i <= length($seq1_no_linefeed)-$K;$i++){
			my $kmer = substr($seq1_no_linefeed,$i,$K);
			
			my $revcom = reverse $kmer;
			$revcom =~ tr/ACGT/TGCA/;
			
			if(exists $jellyfish->{$kmer}){                 # Kmer正向或反向存在于fit_kmer中均加1
				push @kmers_fit,$kmer;
				$count++;
			}elsif(exists $jellyfish->{$revcom}){
				push @kmers_fit,$revcom;
				$count++;
			}
		}
		
		chomp(my $seq2_no_linefeed = $seq2); 			    # 去掉换行符
		for(my $i = 0;$i <= length($seq2_no_linefeed)-$K;$i++){
			my $kmer = substr($seq2_no_linefeed,$i,$K);
			
			my $revcom = reverse $kmer;
			$revcom =~ tr/ACGT/TGCA/;
			
			if(exists $jellyfish->{$kmer}){
				push @kmers_fit,$kmer;                      # 将一对reads的所有uniq kmers都加入至同一数组            
				$count++;
			}elsif(exists $jellyfish->{$revcom}){
				push @kmers_fit,$revcom;
				$count++;
			}
		}

		if($count > ($Length-$K+1)*2*0.8){      # 寻找至少有(125-17+1)*2*0.8=174.4条kmer为uniq kmer的reads进行初始化,2代表一对reads
			for my $kmer (@kmers_fit){          # 将该对reads所有的uniq kmers赋值，进行初始化
				push @{$jellyfish->{$kmer}},$set;
				@{$jellyfish->{$kmer}} = uniq @{$jellyfish->{$kmer}};
			}
			$set++;			
		}
		$count=0;
		@kmers_fit=();
		last LINE if $set > $Group;                         # 找到Group个后即退出循环
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
# 初始化完毕后，重新读取文件
# ===================
my(@kmers,%dataset1,%dataset2,%seq_fq,$fq1,$fq2);
open FL, "<$File" or die "$!";                          # 读取fastq文件列表
while(my $fq = <FL>){
	chomp($fq);
	next if $fq=~/^\s+$/;
	next if $fq=~/(.*?)2\.(fastq|fq)\.gz$/;
	next if $fq=~/(.*?)2\.(fastq|fq)$/;
	
	if($fq=~/(.*?)1\.(fastq|fq)\.gz$/){                 # 打开成对的fastq文件
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
	
	while ( my $header1 = <FQ1> ) {                                      # 打开fq1，读取4行
		my $seq1 = <FQ1>;
		my $qual_header1 = <FQ1>;
		my $qual1 = <FQ1>;
		#die "ERROR: expected \'\@\' but saw $header1" unless substr( $header1, 0, 1 ) eq '@';
		
		my $header2 = <FQ2>;                                             # 打开fq2，读取4行
		my $seq2 = <FQ2>;
		my $qual_header2 = <FQ2>;
		my $qual2 = <FQ2>;
		#die "ERROR: expected \'\@\' but saw $header2" unless substr( $header2, 0, 1 ) eq '@';

		chomp(my $seq1_no_linefeed=$seq1);                               # 去除seq后面的换行符
		$seq_fq{$seq1_no_linefeed}=$header1.$seq1.$qual_header1.$qual1;  # 将fastq文件转化为序列为键，完整的四行fastq文件为值的hash；
		
		chomp(my $seq2_no_linefeed=$seq2);
		$seq_fq{$seq2_no_linefeed}=$header2.$seq2.$qual_header2.$qual2;

		for(my $i=0; $i <= length($seq1_no_linefeed)-$K; $i++){
			my $kmer = substr($seq1_no_linefeed,$i,$K);
			
			my $revcom = reverse $kmer;
			$revcom =~ tr/ACGT/TGCA/;
			
			if(exists $jellyfish->{$kmer}){
				for my $num_set (@{$jellyfish->{$kmer}}){
					$dataset1{$num_set}++;                       # 如果read的某条kmer位于kmer池中，该数据集累计。
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
					$dataset2{$num_set}++;                       # 如果read的某条kmer位于kmer池中，该数据集累计。
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

		my @fit=uniq @fit1,@fit2;                               # 一对reads中，只要有一条满足条件即可

		if($#fit==-1){                                          # 均不满足则输出至remain数据集中
			print REMAIN $seq_fq{$seq1_no_linefeed};
			print REMAIN $seq_fq{$seq2_no_linefeed};
		}
		elsif($#fit==0){                                        # 如果只有一个元素满足，将该read的所有kmers添加至kmer池
			my $file="Group".$fit[0].".fq"; 
			open OUT, ">>$file" or die "$!";
			print OUT $seq_fq{$seq1_no_linefeed};
			print OUT $seq_fq{$seq2_no_linefeed};
			for my $kmer (@kmers){
				push @{$jellyfish->{$kmer}},$fit[0];
				@{$jellyfish->{$kmer}} = uniq @{$jellyfish->{$kmer}};
			}
		}
		else{                                                   # 如果有多个元素满足
			my $first=shift(@fit);                              # 修改kmer池中数据集指向第一个满足条件的数据集
			my $file="Group".$first.".fq";                      # 取第一个满足的数据集的编号
			
			open OUT, ">>$file" or die "$!";
			print OUT $seq_fq{$seq1_no_linefeed};
			print OUT $seq_fq{$seq2_no_linefeed};
			
			unless(system "for i in @fit; do cat Group\$i.fq >> $file; done"){
				for(@fit){
					unlink "Group$_.fq" or warn "$!";           # 运行成功则删除原文件 
					print "Group$_.fq is merged with $file and Group$_.fq has been deleted successfully!\n";
				}
			}

			for my $kmer (@kmers){
				push @{$jellyfish->{$kmer}},$first;
				@{$jellyfish->{$kmer}} = uniq @{$jellyfish->{$kmer}};
			}
			
			for my $num_set (@fit){                             # 判断余下满足条件的值,并去重
				for my $kmer (keys %$jellyfish){
					for(@{$jellyfish->{$kmer}}){
						$_ = $first if $_ == $num_set;
					}
					@{$jellyfish->{$kmer}} = uniq @{$jellyfish->{$kmer}};
				}
			}
		}
		(%seq_fq,%dataset1,%dataset2,@fit,@fit1,@fit2,@kmers)=();                     # 清空以节省内存
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
# 迭代，从Remain.fq中寻找
#==================
my %dataset;
for(my $j=2;$j<=$Iteration;$j++){
	my $p=$j-1;
	open FQ,"<Remain$p.fq" or die "$!";                         # 将Remain1.fq中的reads迭代搜寻
	open REMAIN,">>Remain$j.fq" or die "$!";                    # 将剩下的reads输出至Remain2.fq中
	while ( my $header = <FQ> ) {
		my $seq = <FQ>;
		my $qual_header = <FQ>;
		my $qual = <FQ>;
		#die "ERROR: expected \'\@\' but saw $header" unless substr( $header, 0, 1 ) eq '@';
	
		chomp(my $seq_no_linefeed=$seq);                               # 此处去除seq后面的换行符
		$seq_fq{$seq_no_linefeed}=$header.$seq.$qual_header.$qual;     # 将fastq文件转化为序列为键，完整的四行fastq文件为值的hash；
	
		for(my $i=0;$i<=length($seq_no_linefeed)-$K;$i++){
			my $kmer=substr($seq_no_linefeed,$i,$K);
			
			my $revcom = reverse $kmer;
			$revcom =~ tr/ACGT/TGCA/;
			
			if(exists $jellyfish->{$kmer}){
				for my $num_set (@{$jellyfish->{$kmer}}){
					$dataset{$num_set}++;                       # 如果read的某条kmer位于kmer池中，该数据集累计。
				}
				push @kmers,$kmer;
			}elsif(exists $jellyfish->{$revcom}){
				for my $num_set (@{$jellyfish->{$revcom}}){
					$dataset{$num_set}++;
				}
				push @kmers,$revcom;
			}
		}

		my @fit=grep{$dataset{$_}>$Number}(keys %dataset);      # 如果大于阈值，则将该数据集编号加入@fit中
		
		if($#fit==-1){                                          # 均不满足则输出至remain数据集中
			print REMAIN $seq_fq{$seq_no_linefeed};
		}
		elsif($#fit==0){                                        # 如果只有一个元素满足
			my $file="Group".$fit[0].".fq";                     # 取这个满足的数据集的编号
			open OUT, ">>$file" or die "$!";
			print OUT $seq_fq{$seq_no_linefeed};
			for my $kmer (@kmers){
				push @{$jellyfish->{$kmer}},$fit[0];
				@{$jellyfish->{$kmer}} =  uniq @{$jellyfish->{$kmer}};
			}
		}
		else{                                                   # 如果有多个元素满足
			my $first=shift(@fit);                              # 修改kmer池中数据集指向第一个满足条件的数据集
			my $file="Group".$first.".fq";                      # 取第一个满足的数据集的编号
			
			open OUT, ">>$file" or die "$!";
			print OUT $seq_fq{$seq_no_linefeed};
			
			unless(system "for i in @fit; do cat Group\$i.fq >> $file; done"){
				for(@fit){
					unlink "Group$_.fq" or warn "$!";           # 运行成功则删除原文件 
					print "Group$_.fq is merged with $file and Group$_.fq has been deleted successfully!\n";
				}
			}

			for my $kmer (@kmers){
				push @{$jellyfish->{$kmer}},$first;
				@{$jellyfish->{$kmer}} = uniq @{$jellyfish->{$kmer}};
			}
			
			for my $num_set (@fit){                             # 判断余下满足条件的值,并去重
				for my $kmer (keys %$jellyfish){
					for(@{$jellyfish->{$kmer}}){
						$_ = $first if $_ == $num_set;
					}
					@{$jellyfish->{$kmer}} = uniq @{$jellyfish->{$kmer}};
				}
			}
		}
		(%seq_fq,%dataset,@fit,@kmers)=();                     # 清空以节省内存		
	}
	close FQ;
	close OUT;
	close REMAIN;
	
	unlink "Remain$p.fq" or warn "$!";                         # 删除上次迭代的Remain.fq
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
