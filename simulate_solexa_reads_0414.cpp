#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <map>

using namespace std;

char *input;
int read_length=100;
int coverage=40;
int insertsize_mean=500;
int insertsize_sd=25;
float error_rate=0.01;
float hetersnp_rate=0;
float heterindel_rate=0;
char *output;

const char *MAKE_TIME="2010-03-01";
const char *VERSION="1.0";
const char *AUTHOR="lujianliang";
const char *CONTACT="lujianliang@genomics.org.cn";

void Usage(){
	cout<<"Description:"<<endl;
	cout<<endl<<"\tIt is a program for simulating solexa PE reads,with a series of problems generate by solexa"<<endl;
	cout<<"\tsequencing machine in thought,such as insertsize distribution,error rate,heterozygosis SNP"<<endl;
	cout<<"\tand heterozygosis Indel in diploid. User should set the value of insertsize_mean and"<<endl;
	cout<<"\tinsertsize_sd ,they are the mean value and standard deviation of the normal distribution"<<endl;
	cout<<"\tthat used as the model function when simulating insertsize distribution,usually we set the"<<endl;
	cout<<"\tinsertsize_sd to be 1/20 of the insertsize_mean.The normal distribution function model we"<<endl;
	cout<<"\tused in this program is f(x)=1/考/sqrt(2*pi)/exp((x-米)**2 / (2*考**2)) ,and the insertsize"<<endl;
	cout<<"\tdistribution range is limited in (米-5考ㄛ米+5考), which will cover almost all the data."<<endl;
	cout<<"\tTo simulate illumina error rates on different cycles,we use the function f(x)=0.00001*x**3"<<endl;
	cout<<"\tas a model function,because the error rate of a number of bases at the end is much larger"<<endl;
	cout<<"\tthan other area in a read. User should set the heterozygosis SNP rate and heterozygosis Indel"<<endl;
	cout<<"\trate as also ,if you want to use those function,but remember that heterozygosis SNP rate"<<endl;
	cout<<"\tand heterozygosis Indel rate is only exists in diploid. At last ,you should set another"<<endl;
	cout<<"\tseveral parameters ,read length ,coverage of reads ,input sequence and output prefix,the"<<endl;
	cout<<"\tinput sequence must be set ,because there is not default value."<<endl;
	cout<<endl<<"Program:simulate_solexa_reads"<<endl;
	cout<<"\tCompile Data:\t"<<MAKE_TIME<<endl;
	cout<<"\tAuthor:\t"<<AUTHOR<<endl;
	cout<<"\tVersion:\t"<<VERSION<<endl;
	cout<<"\tContact:\t"<<CONTACT<<endl;
	cout<<endl<<"Usage:\tsimulate_solexa_reads [options]"<<endl;
	cout<<"\t-i 	<string>	input,input reference genome sequence *.fa"<<endl;
	cout<<"\t-l 	<int>		read_len,set read length,read1 and read2 have the same length,default:100"<<endl;
	cout<<"\t-x 	<int>		coverage,set	the sequencing coverage(sometimes called depth),default:40"<<endl;
	cout<<"\t-m 	<int>		insertsize_mean,set the average value of insert size,default:500"<<endl;
	cout<<"\t-v 	<int>		insertsize_sd,set the standard deviation of insert sizes, default:25"<<endl;
	cout<<"\t-e 	<float>		error_rate,set the average error rate over all cycles,default:0.01"<<endl;
	cout<<"\t-s 	<float>		heterSNP_rate,set the heterozygous SNP rate of the diploid genome,default:0"<<endl;
	cout<<"\t-d 	<float>		heterIndel_rate,set the heterozygous indel rate of the diploid genome,default:0"<<endl;
	cout<<"\t-o 	<string>	output,output file prefix default:solexa"<<endl;
	cout<<"\t-h 				help,output help infomation"<<endl;
	cout<<endl<<"Example:"<<endl;
	cout<<endl<<"\t1. perl simulate_solexa_reads.pl -i ref_sequence.fa"<<endl;
	cout<<"\tEvery parameter use the default one."<<endl;
	cout<<"\t2. perl simulate_solexa_reads.pl -i ref_sequence.fa -l 150 -x 20 -o humen_500_100"<<endl;
	cout<<"\tJust set read length and coverage you needed."<<endl;
	cout<<"\t3. perl simulate_solexa_reads.pl -i ref_sequence.fa -o humen -m 600 -v 30 -e 0.01"<<endl;
	cout<<"\tSet insertsize distribution and error rate."<<endl;
	cout<<"\t4. perl simulate_solexa_reads.pl -i ref_sequence.fa -o humen -s 0.001 -d 0.001"<<endl;
	cout<<"\tThe genome is diploid and you want to produce heterozygosis SNPs  heterozygosis Indels in reads."<<endl;
	exit(-1);
}

void Getopt(int argc,char *argv[]){
	int c;
	while ((c=getopt(argc,argv,"i:l:x:m:v:e:s:d:o:h"))!=-1)
	{
		switch(c){
			case 'i': input=optarg;break;
			case 'l': read_length=atoi(optarg);break;
			case 'x': coverage=atoi(optarg);break;
			case 'm': insertsize_mean=atoi(optarg);break;
			case 'v': insertsize_sd=atoi(optarg);break;
			case 'e': error_rate=atof(optarg);break;
			case 's': hetersnp_rate=atof(optarg);break;
			case 'd': heterindel_rate=atof(optarg);break;
			case 'o': output=optarg;break;
			case 'h': Usage();break;
			default: Usage();
		}
	}
}

//simulate the insertsize distribution with the model of normal distribution function
//The insertsize range is limited in (米-5考ㄛ米+5考), which covers almost all the data.
vector <int> insert_distribution(int reads_pair){
	double pi=3.1415926535;
	vector <double> insert;
	vector <int> insert_num;
	double total,temp1;
	int temp2=0,total2=0,num=0;
	for (int i=insertsize_mean-5*insertsize_sd;i<=insertsize_mean+5*insertsize_sd;i++)
	{
		temp1=1/sqrt(2*pi)/insertsize_sd/exp(pow((i-insertsize_mean),2)/(2*pow(insertsize_sd,2)));
		insert.push_back(temp1);
		total+=temp1;
	}
	for (int i=insertsize_mean-5*insertsize_sd;i<=insertsize_mean+5*insertsize_sd;i++){
		temp2=int(0.5+insert[num++]/total*reads_pair);
		insert_num.push_back(temp2);
		total2+=temp2;
	}
	insert_num[5*insertsize_sd]+=reads_pair-total2;
	insert.clear();
	return insert_num;
}

//Simulate illumina error distribution on different cycles,
//with the model function f(x)=0.00001*x**3
vector <double> error_distribution(int rd_pair){
	double basic_error_rate=0.001;
	double total_error=(error_rate-basic_error_rate)*read_length;
	vector <double> error_dist;
	double total;
	for (int i=1;i<=read_length;i++)
	{
		double temp=0.00001*pow(i,3);
		error_dist.push_back(temp);
		total+=temp;
	}
	for (int i=1;i<=read_length;i++)
	{
		double temp=(basic_error_rate+error_dist[i-1]/total*total_error)*rd_pair;
		error_dist[i-1]=temp;
	}
	return error_dist;
}

//Rrealization of error sequencing
char get_match(char base){
	char bases[4][3]={{'T','G','C'},
			 {'A','G','C'},
			 {'A','T','G'},
			 {'A','T','C'}};
//	map<char,int>match;
//	match['A']=0;
//	match['T']=1;
//	match['C']=2;
//	match['G']=3;
	int a;
	char n='N';
	switch (base)
	{
	case 'A': a=0;break;
	case 'T': a=1;break;
	case 'C': a=2;break;
	case 'G': a=3;break;
	default: return n;
	}
	int num=int(rand()%3);
	return bases[a][num];
}

//Produce heterozygous SNPs in multiploid
string Get_snp(string &seq,ofstream &snp,string id){
	int seq_len=seq.size();
	int snp_num=int(seq_len*hetersnp_rate);
	for (int i=0;i<snp_num;i++)
	{
		int index=int(rand()%seq_len);
		snp<<id<<"\t"<<index<<"\t"<<seq[index]<<"\t";
		seq[index]=get_match(seq[index]);
		snp<<seq[index]<<endl;
	}
	return seq;
}

//Getting insert sequence when heterozygous indel rate is bigger than 0
string get_insertion(int num){
	char base[]={'A','T','G','C'};
	string s;
	for (int a=0;a<num;a++)
	{
		int index=int(rand()%4);
		s+=base[index];
	}
	return s;
}

//Produce heterozygous indels in multiploid
string Get_indel(string &seq,ofstream &indel,string id1){
	int seq_len=seq.size();
	int indel_num=int(seq_len*heterindel_rate);
	int array[3]={2,3,6};
	int p=1;
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<indel_num/2/array[i];j++)
		{
			int num=int(rand()%seq_len);
			if (num+p>seq_len)
			{
				j--;
			}else{
				indel<<id1<<"\t"<<"-"<<"\t"<<num+1<<"\t"<<p<<"\t";
				for (int k=0;k<p;k++)
				{
					indel<<seq[num+k];
					seq[num+k]='N';
				}
				indel<<endl;
			}	
		}
		p++;
	}
	p=1;
	vector<int> insert(seq_len);
	string s;
	for (int i=0;i<3;i++)
	{
		for (int j=0;j<indel_num/2/array[i];j++)
		{
			int num=int(rand()%seq_len);
			insert[num]=p;
		}
		p++;
	}
	for (int i=0;i<seq_len;i++)
	{
		if (insert[i]>=1)
		{
			string temp;
			temp=get_insertion(insert[i]);
			s+=temp;
			indel<<id1<<"\t+\t"<<i+1<<"\t"<<insert[i]<<"\t"<<temp<<endl;
			if (seq[i]!='N')
			{
				s+=seq[i];
			}
		}else{
			if (seq[i]!='N')
			{
				s+=seq[i];
			}
		}
	}
	return s;
}

//get reverse_complementary sequence
string reversecomplementary(string read){
	string s;
	for (int i=read.size()-1;i>=0;i--)
	{
		if (read[i]=='A')
		{
			s+='T';
		}else if (read[i]=='T')
		{
			s+='A';
		}else if (read[i]=='G')
		{
			s+='C';
		}else{
			s+='G';
		}
	}
	return s;
}

long long output_reads(vector <int> insert_dist,vector <double> err_dist,string &seq,int seqlen,
	int rd_pair,string id_seq,ofstream &o1,ofstream &o2,ofstream &log3,long long reads_all){
	int reads_count=0;
	int k=0;
	log3<<"Begin to output reads"<<endl;
	srand((unsigned)time(NULL));
	for (int i=insertsize_mean-5*insertsize_sd;i<=insertsize_mean+5*insertsize_sd;i++)
	{	
		if (i<read_length)
		{
			continue;
		}
		if (seqlen<i)
		{
			return reads_count;
		}
		while (insert_dist[k]>0)
		{  
			int pos=int(rand()%seqlen);
			if (pos+i>seqlen)
			{
				continue;
			}
			reads_count++;
			string substr=seq.substr(pos,i);
			string read1=substr.substr(0,read_length);
			string read2=substr.substr(i-read_length,read_length);
			if (insertsize_mean>=1000)
			{
				read1=reversecomplementary(read1);
			}else if (insertsize_mean>0 && insertsize_mean<1000)
			{
				read2=reversecomplementary(read2);
			}else{
				cout<<"Error:insertsize is smaller than 0"<<endl;
				exit(-1);
			}
			int errorpos1[100];
			int errorpos2[100];
			char errorbase1[100];
			char errorbase2[100];
			int error_count1=0;
			int error_count2=0;
			if (error_rate>=0.001)
			{	
				for (int j=0;j<read_length;j++)
				{
					double num=double(rand()%rd_pair);
					if (num<err_dist[j-1])
					{
						errorpos1[error_count1]=j;
						errorbase1[error_count1]=read1[j];
						read1[j]=get_match(read1[j]);
						error_count1++;
					}
				}
				for (int j=0;j<read_length;j++)
				{
					double num=double(rand()%rd_pair);
					if (num<err_dist[j-1])
					{
						errorpos2[error_count2]=j;
						errorbase2[error_count2]=read2[j];
						read2[j]=get_match(read2[j]);
						error_count2++;
					}
				}
			}
			int num1=int(rand()%2);
			if (reads_count%10000==0)
			{
				log3<<"Output "<<reads_count<<" reads"<<endl;
			}
			if (num1==0)
			{
//				o1<<">"<<id_seq<<"_"<<insertsize_mean<<"_"<<read_length<<"_"<<pos<<"_"<<reads_count<<"_"<<1<<endl<<read1<<endl;
				o1<<">read_"<<insertsize_mean<<"_"<<reads_count+reads_all<<"_"<<1<<" "<<id_seq<<" "<<pos+1<<" "<<read_length<<" ";
				if (error_count1>0)
				{
					for (int a=0;a<error_count1;a++)
					{
						o1<<errorpos1[a]+1<<","<<errorbase1[a]<<";";
					}
				}
				o1<<endl<<read1<<endl;
//				o2<<">"<<id_seq<<"_"<<insertsize_mean<<"_"<<read_length<<"_"<<pos<<"_"<<reads_count<<"_"<<2<<endl<<read2<<endl;
				o2<<">read_"<<insertsize_mean<<"_"<<reads_count+reads_all<<"_"<<2<<" "<<id_seq<<" "<<pos+i-read_length+1<<" "<<read_length<<" ";
				if (error_count2>0)
				{
					for (int a=0;a<error_count2;a++)
					{
						o2<<errorpos2[a]+1<<","<<errorbase2[a]<<";";
					}
				}
				o2<<endl<<read2<<endl;
			}else{
//				o1<<">"<<id_seq<<"_"<<insertsize_mean<<"_"<<read_length<<"_"<<pos<<"_"<<reads_count<<"_"<<1<<endl<<read2<<endl;
				o1<<">read_"<<insertsize_mean<<"_"<<reads_count+reads_all<<"_"<<1<<" "<<id_seq<<" "<<pos+i-read_length+1<<" "<<read_length<<" ";
				if (error_count2>0)
				{
					for (int a=0;a<error_count2;a++)
					{
						o1<<errorpos2[a]+1<<","<<errorbase2[a]<<";";
					}
				}
				o1<<endl<<read2<<endl;
//				o2<<">"<<id_seq<<"_"<<insertsize_mean<<"_"<<read_length<<"_"<<pos<<"_"<<reads_count<<"_"<<2<<endl<<read1<<endl;
				o2<<">read_"<<insertsize_mean<<"_"<<reads_count+reads_all<<"_"<<2<<" "<<id_seq<<" "<<pos+1<<" "<<read_length<<" ";
				if (error_count1>0)
				{
					for (int a=0;a<error_count1;a++)
					{
						o2<<errorpos1[a]+1<<","<<errorbase1[a]<<";";
					}
				}
				o2<<endl<<read1<<endl;
			}
			insert_dist[k]--;
		}
		k++;
	}
	log3<<"Finish output reads"<<endl;
	return reads_count;
}

long long get_reads(string id,string &sequ,ofstream &of1,ofstream &of2,ofstream &log2,long long read_genome,ofstream &snp,ofstream &indel){
	string sequence;
	long long readonchr=0;
	if (sequ.size()<insertsize_mean)
	{
		return 0;
	}
	for (int i=0;i<sequ.size();i++)
	{
		switch (sequ[i])
		{
			case 'a': sequence+='A';break;
			case 't': sequence+='T';break;
			case 'c': sequence+='C';break;
			case 'g': sequence+='G';break;
			case 'A': sequence+='A';break;
			case 'T': sequence+='T';break;
			case 'G': sequence+='G';break;
			case 'C': sequence+='C';break;
			default : break;
		}
	}
	int sequence_length=sequence.size();
	long long reads_pair_num=0;

	if ((hetersnp_rate>0) || (heterindel_rate>0))
	{
		reads_pair_num=(long long)sequence_length*coverage/read_length/2/2;
	}else{
		reads_pair_num=(long long)sequence_length*coverage/read_length/2;
	}
	vector<int>insertsize_distribution;
	log2<<"Begin to simulate insertsize distrubution"<<endl;
	insertsize_distribution=insert_distribution(reads_pair_num);
	log2<<"Finish simulating insertsize distrubution"<<endl;
	vector <double> circle_error_rate;
	if (error_rate>0)
	{
		if (error_rate<0.001)
		{
			cout<<"Error:error_rate is smaller than basic error rate 0.001,you'd better set the value" 
			<<"either bigger or 0"<<endl;
			exit(-1);
		}else{
			log2<<"Begin to simulate error distrubution"<<endl;
			circle_error_rate=error_distribution(reads_pair_num);
			log2<<"Finish to simulate error distrubution"<<endl;
		}
	}
	readonchr=output_reads(insertsize_distribution,circle_error_rate,sequence,sequence_length,reads_pair_num,
		id,of1,of2,log2,read_genome);
	if (hetersnp_rate>0 || heterindel_rate>0) //heterozygous SNP and heterozygous indel exists is diploid
	{
		if (hetersnp_rate>0)
		{
			sequence=Get_snp(sequence,snp,id);
		}
		if (heterindel_rate>0)
		{
			log2<<"Begin to simulate indel"<<endl;
			sequence=Get_indel(sequence,indel,id);
			log2<<"Finish to simulate indel"<<endl;
		}
		sequence_length=sequence.size();
		reads_pair_num=(long long)sequence_length*coverage/read_length/2/2;
		vector <int> insertsize_distribution;
		log2<<"Begin to simulate insertsize distrubution"<<endl;
		insertsize_distribution=insert_distribution(reads_pair_num);
		log2<<"Finish simulating insertsize distrubution"<<endl;
		vector <double> circle_error_rate;
		if (error_rate>0)
		{
			log2<<"Begin to simulate error distrubution"<<endl;
			circle_error_rate=error_distribution(reads_pair_num);
			log2<<"Finish to simulate error distrubution"<<endl;
		}
		long long readonchr2=read_genome+readonchr;
		readonchr+=output_reads(insertsize_distribution,circle_error_rate,sequence,sequence_length,reads_pair_num,
			id,of1,of2,log2,readonchr2);
	}
	return readonchr;
}

void Get_genome(ifstream &inf,ofstream &outf1,ofstream &outf2,ofstream &log1,ofstream &snp,ofstream &indel){
	string line,id,seq;
	long long readINgenome=0;
	while (getline(inf,line,'\n'))
	{
		if (line[0]=='>')
		{
			if (seq!="")
			{	
				log1<<"Have finished reading scaffold "<<id<<endl;
				readINgenome+=get_reads(id,seq,outf1,outf2,log1,readINgenome,snp,indel);
				seq="";
			}
			line.erase(0,1);
//			id=line;
			int pos=line.find(" ");
			line=line.substr(0,pos);
			id=line;
		}else{
			seq+=line;
		}		
	}
	log1<<"Have finished reading scaffold "<<id<<endl;
	readINgenome+=get_reads(id,seq,outf1,outf2,log1,readINgenome,snp,indel);
}

int main(int argc, char *argv[])
{
	if (argc==1)
	{
		Usage();
	}
	Getopt(argc,argv);
	ifstream infile;
	ofstream outfile1;
	ofstream outfile2;
	ofstream log;
	infile.open(input);
	ofstream snp,indel;
	if (!infile)
	{
		cerr<<"Error:unable to open input file:"<<input<<endl;
		exit(-1);
	}
	
	if (!output)
	{
		char out1[100];
		char out2[100];
		char out3[100];
		sprintf(out1,"%s%d%s%d%s","solexa_",read_length,"_",insertsize_mean,"_1.fa");
		sprintf(out2,"%s%d%s%d%s","solexa_",read_length,"_",insertsize_mean,"_2.fa");
		sprintf(out3,"%s%d%s%d%s","solexa_",read_length,"_",insertsize_mean,".log");
		outfile1.open(out1);
		outfile2.open(out2);
		log.open(out3);
		if (hetersnp_rate>0)
	         {
//        	         char snp_array[]="snp_solexa.lis";
			 char snp_array[100];
                	 sprintf(snp_array,"%s%d%s%d%s","snp_solexa_",read_length,"_",insertsize_mean,".lis");
                 	 snp.open(snp_array);
        	 }
        	 if (heterindel_rate>0)
        	 {
                	 char indel_array[]="indel_solexa.lis";
                 	sprintf(indel_array,"%s%d%s%d%s","indel_solexa_",read_length,"_",insertsize_mean,".lis");
                	 indel.open(indel_array);
         	}

	}else{
		char out1[500];char out2[500];char out3[500];
/*		string s1(output);
		string s2(output);
		string s3(output);
		s1=s1+"_"+read_length+"_"+insertsize_mean+"_1.fa";
		s2=s2+"_"+read_length+"_"+insertsize_mean+"_2.fa";
		s3=s3+"_"+read_length+"_"+insertsize_mean+".log";
		strcpy(out1,s1.c_str());
		strcpy(out2,s2.c_str());
		strcpy(out3,s3.c_str());
*/		
		sprintf(out1,"%s%s%d%s%d%s",output,"_",read_length,"_",insertsize_mean,"_1.fa");
		sprintf(out2,"%s%s%d%s%d%s",output,"_",read_length,"_",insertsize_mean,"_2.fa");
		sprintf(out3,"%s%s%d%s%d%s",output,"_",read_length,"_",insertsize_mean,".log");
		outfile1.open(out1);
		outfile2.open(out2);
		log.open(out3);
		if (hetersnp_rate>0)
	         {
        	         char snp_array[100];
               		 sprintf(snp_array,"%s%s%s%d%s%d%s","snp_",output,"_",read_length,"_",insertsize_mean,".lis");
               		 snp.open(snp_array);
        	 }
         	if (heterindel_rate>0)
        	 {
                	 char indel_array[100];
                	 sprintf(indel_array,"%s%s%s%d%s%d%s","indel_",output,"_",read_length,"_",insertsize_mean,".lis");
                	 indel.open(indel_array);
        	 }

	}
//float hetersnp_rate=0;
//float heterindel_rate=0;
/*	ofstream snp,indel;
	if (hetersnp_rate>0)
	{
		char snp_array[100];
		sprintf(snp_array,"%s%d%s","snp_",insertsize_mean,".lis");
		snp.open(snp_array);
	}
	if (heterindel_rate>0)
	{
		char indel_array[100];
		sprintf(indel_array,"%s%d%s","indel_",insertsize_mean,".lis");
		indel.open(indel_array);
	}
*/
	if (!outfile1 || !outfile2 || !log)
	{
		if (!outfile1 || !outfile2)
		{
			cerr<<"Error:unable to open output file."<<endl;
			exit(1);
		}else{
			cerr<<"Error:unable to open output *log file."<<endl;
			exit(1);
		}
	}
	Get_genome(infile,outfile1,outfile2,log,snp,indel);
	return 0;
}
