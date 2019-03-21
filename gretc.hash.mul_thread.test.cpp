#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string.h>
#include <bitset>
#include <vector>
#include "gzstream.h"
#include <algorithm>
#include <pthread.h>

using namespace std;

#define MAX_GROUP_NUM 10
#define DEBUG true
/*
#define KMER_SIZE 17
#define KMER_BI_SIZE 34
#define HASH_SIZE 15
#define HASH_BI_SIZE 30
#define REST_KMER_SIZE 2
#define REST_KMER_BI_SIZE 4
#define MAX_READS_LENGTH 500
#define MAX_READS_BI_LENGTH 1000
#define BUCKET_SIZE 1073741824
*/


#define KMER_SIZE 17
#define KMER_BI_SIZE 34
#define HASH_SIZE 12
#define HASH_BI_SIZE 24
#define REST_KMER_SIZE 5
#define REST_KMER_BI_SIZE 10
#define MAX_READS_LENGTH 500
#define MAX_READS_BI_LENGTH 1000
#define BUCKET_SIZE 16777216


typedef struct GroupList
{
	unsigned int group_id;
	struct GroupList *next;
	GroupList(unsigned int id);
}GList;

GroupList::GroupList(unsigned int id)
{
	group_id = id;
	next = NULL;
}

//四叉树节点
typedef struct kmer
{
	struct kmer *base[4];						//指针数组，分别指向actg,base[0]=*a,base[1]=*c,base[2]=*t,base[3]=*g
	struct GroupList *group;
	bitset<REST_KMER_BI_SIZE> k;
	kmer();
}Kmer;

kmer::kmer()								//结构体的默认值
{
	memset(base,0,sizeof(base));
	group = NULL;
}


//全局变量
bitset<MAX_GROUP_NUM> group_infor;						//用于标记每个group的状态，1表示已启用，0表示未启用
bitset<REST_KMER_BI_SIZE> result;						//用于存储两个bitset异或的结果
bitset<REST_KMER_BI_SIZE> kmer_complement;					//用于kmer求互补链
bitset<MAX_READS_BI_LENGTH> reads_complement;					//用于reads求互补链
bitset<MAX_GROUP_NUM> flush_group_bitset;					//用于更新group的状态
float overlap_threshold = 0.2;
float initialize_threshold = 0.8;
string output_path = "./";							//输出路径
Kmer ** hash = NULL;								//存储uniq-kmer的后缀树
vector<string> current_file;							//分reads过程，当前需要分的文件
igzstream reads_file1;								//多线程中的reads1文件
igzstream reads_file2;								//多线程中的reads2文件
long finish_reads_num = 0;							//记录当前完成了多少reads
long groups_base_num[MAX_GROUP_NUM+1];						//记录每个group的碱基数
size_t hash_idx = 0;								//hash数组的小标
vector<size_t> final_group_st;							//final_group的size_t形式
vector<unsigned int> final_group_ui;						//final_group的unsigned形式

int calculate_thread_num = 1;
int flush_thread_num = 1;

pthread_mutex_t reads_lock;							//多线程读入文件锁，锁住reads_file1,reads_file2
pthread_mutex_t assign_lock;							//多线程统计输出的锁
pthread_mutex_t flush_lock;							//多线程刷新锁,锁住hash_idx

size_t group_num_id = 100;							//GROUP的数量，size_t是unsigned int类型，作为数据下标
unsigned int group_num_no = 100;						//作为group的最大编号

//将KMER序列转化成二进制，每个字节用两个字符表示：A 00;C 01;T 10;G 11;
void kmer_to_bitset(const string &str,bitset<REST_KMER_BI_SIZE> &kmer_bit,bitset<HASH_BI_SIZE> &hash_bit)
{
	int kmer_i = REST_KMER_BI_SIZE - 1;
	int hash_i = HASH_BI_SIZE - 1;
	int count = 0;
	for(string::size_type index = 0;index != str.size();++index)
	{
		bitset<8> char_bit(str[index]);
		if(count < HASH_SIZE)
		{
			hash_bit[hash_i--] = char_bit[2];					//字符二进制的第6位
			hash_bit[hash_i--] = char_bit[1];					//字符二进制的第7位
		}else
		{
			kmer_bit[kmer_i--] = char_bit[2];
			kmer_bit[kmer_i--] = char_bit[1];
		}
		++count;
	}
}

//将KMER序列转化成互补链，并转成二进制，每个字节用两个字符表示：A 00;C 01;T 10;G 11;
/*void kmer_to_complementary_bitset(const string &str , bitset<KMER_BI_SIZE> &kmer_bit,bitset<HASH_BI_SIZE> &hash_bit)
{
	unsigned int kmer_i = 0;
	unsigned int hash_i = 0;
	for(string::size_type index = str.size();index !=0 ;){
		--index;
		if(index)
		bitset<8> char_bit(str[index]);
		bit[i++] = char_bit[2];					//字符二进制的第6位
		bit[i++] = char_bit[1];					//字符二进制的第7位
	}
	bit = bit ^ kmer_complement;
}*/

//将reads序列转化成二进制，每个字节用两个字符表示：A 00;C 01;T 10;G 11;
void reads_to_bitset(const string &str , bitset<MAX_READS_BI_LENGTH> &bit)
{
	int i = MAX_READS_BI_LENGTH - 1;
	for(string::size_type index = 0;index != str.size();++index){
		bitset<8> char_bit(str[index]);
		bit[i--] = char_bit[2];					//字符二进制的第6位
		bit[i--] = char_bit[1];					//字符二进制的第7位
	}
}

//将序列转化成互补链，并转成二进制，每个字节用两个字符表示：A 00;C 01;T 10;G 11;
void reads_to_complementary_bitset(const string &str , bitset<MAX_READS_BI_LENGTH> &bit)
{
	int i = MAX_READS_BI_LENGTH - 1;
	for(string::size_type index = str.size();index !=0 ;){
		--index;
		bitset<8> char_bit(str[index]);
		bit[i--] = char_bit[2];					//字符二进制的第6位
		bit[i--] = char_bit[1];					//字符二进制的第7位
	}
	//cout << bit << endl;
	//cout << reads_complement << endl;
	bit = bit ^ reads_complement;
	//cout << bit << endl;
}

//测试用，输出bitset的内容；从0 - bitset.length
void print_bitset(bitset<REST_KMER_BI_SIZE> &bi)
{
	int index = 0;
	for(index = 0; index != REST_KMER_BI_SIZE;++index)
	{
		cout << bi[index];
	}
	cout << endl;
}

//返回下一个base的指针对应的数组下标
void find_next_index(int position,const bitset<REST_KMER_BI_SIZE> &bi,int &index){			//position从0开始
	index = 0;
	if(bi[2 * position])							//高位
	{
		index += 2;
	}

	if(bi[2 * position + 1])
	{
		++index;							//低位
	}
	//cout << index <<endl;
}

//递归查找，插入节点
void insert_to(const bitset<REST_KMER_BI_SIZE> &bi,Kmer *new_kmer,Kmer *temp_kmer,int &position,int &index)
{
	find_next_index(position,bi,index);
	Kmer *base_kmer = temp_kmer->base[index];

	if(base_kmer)//base不为空
	{
		temp_kmer = base_kmer;
		++position;
		insert_to(bi,new_kmer,temp_kmer,position,index);
	}else
	{
		temp_kmer->base[index] = new_kmer;
	}
}


//刷新一个节点，函数输入是全局变量 vectot<size_t> final_group_st，保持链表的有序
//两个链表整合。保留final_group_st的第一个group，去除其余group
void flush_tree_node(Kmer *node)
{
	//1. 先插入final_group_st 的首节点,如果已经包含，则不必处理
	bool flag = true;//true说明已经包含该group
	GList * temp_group = node->group;
	GList *last_group = NULL;

	while(temp_group)
	{
		if(temp_group->group_id < final_group_ui[0])
		{
			last_group = temp_group;
			temp_group = temp_group->next;
		}
		else if(temp_group->group_id == final_group_ui[0])
		{
			flag = false;
			break;	
		}else
		{
			break;
		}
	}
	
	//在链表头部插入;一种是尚未有链表，另一种是group_id 小于表头
	if(flag)
	{
		GList *new_group = new GList(final_group_ui[0]);
		new_group->next = temp_group;
		if(last_group)//在中部或尾部插入
		{
			last_group->next = new_group;
		}else//在头部插入
		{
			node->group = new_group;
		}
	}

	//2. 删掉后面的节点
	for(vector<size_t>::size_type id = 1;id != final_group_ui.size();++id)
	{
		while(temp_group && temp_group->group_id < final_group_ui[id])
		{
			last_group = temp_group;
			temp_group = temp_group->next;
		}

		if(temp_group)
		{
			//不存在删除表头的情况
			if(temp_group->group_id == final_group_ui[id])//删掉该节点
			{
				last_group->next = temp_group->next;
				temp_group->next = NULL;
				delete temp_group;
				temp_group = last_group->next;
			}else
			{
				continue;				
			}
		}
		else//说明已到group list的尾部
		{
			break;
		}
	}
}

//递归深度优先遍历四叉树,如果flag=1.刷新group(根据全局变量bitset<MAX_GROUP_NUM> flush_group_bitset)
void tranvers_tree(const Kmer *head)
{
	if(head->base[0])
	{
		flush_tree_node(head->base[0]);
//		cout << head->base[0]->k << "+" << head->base[0]->group << endl;
		tranvers_tree(head->base[0]);
	}
	if(head->base[1])
	{
		flush_tree_node(head->base[1]);
//		cout << head->base[1]->k << "+" << head->base[1]->group << endl;
		tranvers_tree(head->base[1]);
	}
	if(head->base[2])
	{
		flush_tree_node(head->base[2]);
//		cout << head->base[2]->k << "+" << head->base[2]->group << endl;
		tranvers_tree(head->base[2]);
	}
	if(head->base[3])
	{
		flush_tree_node(head->base[3]);
//		cout << head->base[3]->k << "+" << head->base[3]->group << endl;
		tranvers_tree(head->base[3]);
	}
	return;
}

//遍历整个散列表,调用tranvers_tree
void *tranvers_hash(void *ar)
{
	size_t temp_idx;
	while(true)
	{
		pthread_mutex_lock(&flush_lock);
		temp_idx = hash_idx++;
		pthread_mutex_unlock(&flush_lock);

		//退出线程
		if(temp_idx >= BUCKET_SIZE){break;}

		if(hash[temp_idx])
		{
			//cout << ix << ": " << endl;
			tranvers_tree(hash[temp_idx]);
		}
	}
}

//根据给定的分割符切割字符串，并返回数组。且字符串的长度小于500bp,大于500bp的分割成多条
void split_string(string str,vector<string> &vecstr,string deli,bool split_to_max_length)
{
        string::size_type position=0;
        while((position=str.find(deli)) != string::npos)
        {
                if(position == 0)
                {
                        str = str.substr(position+1,str.size());
                        continue;
                }else{
                        vecstr.push_back(str.substr(0,position));
                        str = str.substr(position + 1);				//到字符结尾
                }
        }
        if(str.size())vecstr.push_back(str);
	if(!split_to_max_length)return;						//如果不需要切割，则返回

	//对大于500bp的进行切割
	for(vector<string>::size_type idx =0; idx != vecstr.size();++idx)
	{
		if(vecstr[idx].size() > MAX_READS_LENGTH)
		{
			vecstr.push_back(vecstr[idx].substr(MAX_READS_LENGTH - KMER_SIZE + 1));	
			vecstr[idx] = vecstr[idx].substr(0,MAX_READS_LENGTH);
		}
	}

}

//判断两个bitset是否相同,相同返回ture，不同返回false
bool is_equal(const bitset<REST_KMER_BI_SIZE> &b1, const bitset<REST_KMER_BI_SIZE> &b2)
{
	result = b1 ^ b2;	//异或，不一样时为1
	return result.none();	//有1说明不相同，none()没有1
}

//通过kmer的二进制序列，查找在四叉树中的节点
Kmer * seek_in_tree(const bitset<REST_KMER_BI_SIZE> &bi,const Kmer *head_kmer,int &position,int &index)
{
	find_next_index(position,bi,index);
	Kmer *next_kmer = head_kmer->base[index];	//next_kmer是栈内存,可用时间换取
	if(!next_kmer) return NULL;
	if(is_equal(bi,next_kmer->k))
	{
		return next_kmer;
	}else
	{
		return seek_in_tree(bi,next_kmer,++position,index);
	}
}

//统计每条reads的uniq_kmer数，并保存uniq———kmer的地址。flags = ture表示正义链。flags = false 代表互补
int stat_reads_uniq_kmer(const string &reads,vector<Kmer *> * &uniq_kmer_address,bool flag)
{
	int position;					//用于seek_in_tree 引用
	int index;					//用于seek_in_tree 引用
	string deli = "N";				//sequence分隔符
	vector<string> fragment;			//subreads
	fragment.reserve(8);				//two reads 可以设置小一些
	bitset<MAX_READS_BI_LENGTH> reads_bit;		//reads 二进制,所有reads都转化到该变量
	string kmer_seq;				//kmer 序列

	//计算reads的uniq_kmer数
	split_string(reads,fragment,deli,true);		//将read1的N处打断
	//for(vector<string>::size_type v_index = 0;v_index != fragment.size();v_index++){cout << fragment[v_index] << endl;}
	for(vector<string>::size_type v_index = 0;v_index != fragment.size();v_index++)
	{
		string::size_type fragment_size = fragment[v_index].size();
		if(fragment_size < KMER_SIZE)continue;
	
		unsigned int kmer_num = fragment_size - KMER_SIZE + 1;							//kmer数量
		//reads 转化成二进制
		if(flag)
		{
			reads_to_bitset(fragment[v_index],reads_bit);
		}else
		{
			reads_to_complementary_bitset(fragment[v_index],reads_bit);
		}
		bitset<HASH_BI_SIZE> hash_bi(reads_bit.to_string(),0,HASH_BI_SIZE);					//REST_KMER 二进制序列
		bitset<REST_KMER_BI_SIZE> kmer_bi(reads_bit.to_string(),HASH_BI_SIZE,REST_KMER_BI_SIZE);		//HASH 二进制序列
		int hash_point = MAX_READS_BI_LENGTH - HASH_BI_SIZE - 1;
		int kmer_point = MAX_READS_BI_LENGTH - KMER_BI_SIZE - 1;

		while(kmer_num-- != 1)
		{
			unsigned long h_index = hash_bi.to_ulong();
			if(hash[h_index])
			{
				position = 0;
				Kmer * temp_node = seek_in_tree(kmer_bi,hash[h_index],position,index);
				if(temp_node){uniq_kmer_address->push_back(temp_node);}
			}
			hash_bi <<= 2;
			hash_bi[1] = reads_bit[hash_point--];
			hash_bi[0] = reads_bit[hash_point--];
			kmer_bi <<= 2;
			kmer_bi[1] = reads_bit[kmer_point--];
			kmer_bi[0] = reads_bit[kmer_point--];
		}
	}
	fragment.clear();
	return (int)uniq_kmer_address->size();

}

//统计两条reads的uniq_kmer条数，并保存uniq_kmer的地址.因为可容许的reads长度长，用于初始化group信息
//bug:这里需要注意，设定的PEreads是短片段的
int stat_pe_uniq_kmer(const string &reads1,const string &reads2,vector<Kmer *> * &uniq_kmer_address)
{
	int uniq_kmer_num = 0;
	uniq_kmer_address->clear();			//使用前先清空
	uniq_kmer_num = stat_reads_uniq_kmer(reads1,uniq_kmer_address,true) + stat_reads_uniq_kmer(reads2,uniq_kmer_address,false);
	return uniq_kmer_num;
}

//小写转大写
char op(char ch)
{
	if(ch>='a'&&ch<='z')
		return ch-32;
	else
		return ch;
}

//先把小写转化成大写，去掉seq的空白字符，换行字符。只留下ATCGN
void string_clear(const string &str,string &new_str)
{
	string temp_str;
	temp_str.resize(str.size());
	transform(str.begin(),str.end(),temp_str.begin(),op);

	string bases = "ATCGN";
	string::size_type begin_loc = 0;
	string::size_type end_loc = 0;
	while((end_loc = temp_str.find_first_not_of(bases,begin_loc)) != string::npos)
	{
		new_str.append(temp_str.substr(begin_loc,end_loc - begin_loc));
		begin_loc = end_loc + 1;
	}
	new_str.append(temp_str.substr(begin_loc));
	
}

//初始化group信息
void initialize_group(vector<string> &file_infor)
{
	int reads_type = atoi(file_infor[0].c_str());			//1为SE，2为PE
	string library_name = file_infor[1];				//文库编号。相同文库会合并
	string file_format = file_infor[2];				//文件类型，fasta/fastq
	string reads1_file_name = file_infor[3];			//reads1文件名
	string reads2_file_name = "" ;					//reads2文件名

	vector<Kmer *> *uniq_kmer_address = new vector<Kmer *>;		//为兼容大的contigs初始化，使用堆内存
	uniq_kmer_address->reserve(20000);				//uniq_kmer_address 设定为最多的情况
	if(reads_type == 2)						//PE reads
	{
		reads2_file_name = file_infor[4];
		igzstream reads1_i_file,reads2_i_file;
		reads1_i_file.open(reads1_file_name.c_str());
		reads2_i_file.open(reads2_file_name.c_str());
		if(!reads1_i_file){cerr << "fail to open input file : " << reads1_i_file << endl;exit(1);}
		if(!reads2_i_file){cerr << "fail to open input file : " << reads2_i_file << endl;exit(1);}
		
		unsigned int initialized_group = 0;				//已经初始化的group个数
		int uniq_kmer_num = 0;					//uniq_kmer 的个数
		while(initialized_group != group_num_no && !reads1_i_file.eof())
		{
			string reads1_id,reads1_seq,reads1_plus,reads1_qual;
			string reads2_id,reads2_seq,reads2_plus,reads2_qual;
	
			if(!file_format.compare("fastq"))
			{
				getline(reads1_i_file,reads1_id);
				getline(reads1_i_file,reads1_seq);
				getline(reads1_i_file,reads1_plus);
				getline(reads1_i_file,reads1_qual);

				getline(reads2_i_file,reads2_id);
				getline(reads2_i_file,reads2_seq);
				getline(reads2_i_file,reads2_plus);
				getline(reads2_i_file,reads2_qual);
			}else if(!file_format.compare("fasta"))
			{
				getline(reads1_i_file,reads1_id);
				getline(reads1_i_file,reads1_seq,'>');
				getline(reads2_i_file,reads2_id);
				getline(reads2_i_file,reads2_seq,'>');
			}else
			{
				cerr << "reads file format must fastq or fasta." <<endl;
				exit(0);
			}

			//将序列转化成正确的格式，除ATCGN以外的字符都会删除
			string new_reads1_seq;
			string new_reads2_seq;
			new_reads1_seq.reserve(reads1_seq.size());
			new_reads2_seq.reserve(reads2_seq.size());

			string_clear(reads1_seq,new_reads1_seq);
			string_clear(reads2_seq,new_reads2_seq);
			//cout << reads1_seq << endl;
			//cout << new_reads1_seq << endl;
			//cout << reads2_seq << endl;
			//cout << new_reads2_seq << endl;

			int reads_kmer_num = (reads1_seq.length() - KMER_SIZE + 1) * 2;
			uniq_kmer_num = stat_pe_uniq_kmer(new_reads1_seq,new_reads2_seq,uniq_kmer_address);//reads2 取互补
			if(uniq_kmer_num < initialize_threshold * reads_kmer_num)
			{
				uniq_kmer_num = stat_pe_uniq_kmer(new_reads2_seq,new_reads1_seq,uniq_kmer_address);//reads1 取互补链
			}

			if(uniq_kmer_num >= 1)//该reads的uniq kmer数大于设定的阈值
			//if(uniq_kmer_num >= initialize_threshold * reads_kmer_num)//该reads的uniq kmer数大于设定的阈值
			{
				for(vector<Kmer *>::size_type idx = 0;idx != uniq_kmer_address->size();++idx)
				{
					GList *grp_head = uniq_kmer_address->at(idx)->group;
					GList *grp = new GList(initialized_group);
					if(!grp_head)
					{
						uniq_kmer_address->at(idx)->group = grp;
						continue;
					}
					while(grp_head->next){grp_head = grp_head->next;}//因为group是从小到大，所以是有序的
					grp_head->next = grp;
				}
				if(DEBUG)cout << initialized_group << endl;
				group_infor[initialized_group] =  true;
				++initialized_group;
				
			}
		}
		reads1_i_file.close();
		reads2_i_file.close();
		delete uniq_kmer_address;
		uniq_kmer_address = NULL;
		if(initialized_group != group_num_no){cerr << "group not enouge." << endl;exit(0);}

	}else if(reads_type == 1)					//SE reads
	{
		igzstream reads1_i_file;
		reads1_i_file.open(reads1_file_name.c_str());
		if(!reads1_i_file){cerr << "fail to open input file : " << reads1_i_file << endl;exit(1);}
		
		int initialized_group = 0;				//已经初始化的group个数
		int uniq_kmer_num = 0;					//uniq_kmer 的个数
		while(initialized_group != group_num_no && !reads1_i_file.eof())
		{
			string reads1_id,reads1_seq,reads1_plus,reads1_qual;
	
			if(!file_format.compare("fastq"))
			{
				getline(reads1_i_file,reads1_id);
				getline(reads1_i_file,reads1_seq);
				getline(reads1_i_file,reads1_plus);
				getline(reads1_i_file,reads1_qual);
			}else if(!file_format.compare("fasta"))
			{
				getline(reads1_i_file,reads1_id);
				getline(reads1_i_file,reads1_seq,'>');
			}else
			{
				cerr << "reads file format must fastq or fasta." <<endl;
				exit(0);
			}

			//将序列转化成正确的格式，除ATCGN以外的字符都会删除
			string new_reads1_seq;
			new_reads1_seq.reserve(reads1_seq.size());
			string_clear(reads1_seq,new_reads1_seq);
			int reads_kmer_num = reads1_seq.length() - KMER_SIZE + 1;
			uniq_kmer_address->clear();
			uniq_kmer_num = stat_reads_uniq_kmer(new_reads1_seq,uniq_kmer_address,true);//reads1 正链
			if(uniq_kmer_num < initialize_threshold * reads_kmer_num)
			{
				uniq_kmer_address->clear();
				uniq_kmer_num = stat_reads_uniq_kmer(new_reads1_seq,uniq_kmer_address,false);//reads1 取互补链
			}

			//cout << uniq_kmer_num << "\t" << reads_kmer_num << "\t" << (float)uniq_kmer_num/reads_kmer_num<< endl;
			if(uniq_kmer_num >= 1)//该reads的uniq kmer数大于设定的阈值
			//if(uniq_kmer_num >= initialize_threshold * reads_kmer_num)//该reads的uniq kmer数大于设定的阈值
			{
				for(vector<Kmer *>::size_type idx = 0;idx != uniq_kmer_address->size();++idx)
				{
					GList *grp_head = uniq_kmer_address->at(idx)->group;
					GList *grp = new GList(initialized_group);
					if(!grp_head)
					{
						uniq_kmer_address->at(idx)->group = grp;
						continue;
					}
					while(grp_head->next){grp_head = grp_head->next;}//因为group是从小到大，所以是有序的
					grp_head->next = grp;
			 	}
				group_infor[initialized_group] =  true;
				++initialized_group;
				if(DEBUG)cout << initialized_group << endl;
				
			}
		}
		reads1_i_file.close();
		delete uniq_kmer_address;
		uniq_kmer_address = NULL;
		if(initialized_group != group_num_no){cerr << "group not enouge." << endl;exit(0);}
		
	}else								//error reads
	{
		cerr << "ERROR: reads type error.must be 1 or 2."<<endl;
		exit(0);
	}

}

//文件复制 file2 append to file1,and delete file2
void copyfile(string &f1_name,string &f2_name)
{
	if(DEBUG)cout << "merge " << endl << f1_name << endl <<  f2_name << endl;
        igzstream file2;
        file2.open(f2_name.c_str());
        ogzstream file1(f1_name.c_str(),std::ios::app);
        file1 << file2.rdbuf();
        file2.close();
        file1.close();
        remove(f2_name.c_str());
}

//文件合并
void merge_file(int reads_type,string &lib_name,string &file_format)
{
	char temp_group[5];
	string group_no;

	if(reads_type == 2)
	{
		
		sprintf(temp_group,"%d",final_group_st[0]);
		group_no = (string)temp_group;

		string file1_name;
		file1_name.reserve(output_path.size() + 50);
		file1_name.append(output_path);
		file1_name.append("/");
		file1_name.append(lib_name);
		file1_name.append(".group_");
		file1_name.append(group_no);
		file1_name.append(".1.");
		file1_name.append(file_format);
		file1_name.append(".gz");
				
		string file2_name;
		file2_name.reserve(output_path.size() + 50);
		file2_name.append(output_path);
		file2_name.append("/");
		file2_name.append(lib_name);
		file2_name.append(".group_");
		file2_name.append(group_no);
		file2_name.append(".2.");
		file2_name.append(file_format);
		file2_name.append(".gz");

		for(vector<size_t>::size_type id = 1;id != final_group_st.size();++id)
		{
			sprintf(temp_group,"%d",final_group_st[id]);
			group_no = (string)temp_group;

			string file3_name;
			file3_name.reserve(output_path.size() + 50);
			file3_name.append(output_path);
			file3_name.append("/");
			file3_name.append(lib_name);
			file3_name.append(".group_");
			file3_name.append(group_no);
			file3_name.append(".1.");
			file3_name.append(file_format);
			file3_name.append(".gz");
				
			string file4_name;
			file4_name.reserve(output_path.size() + 50);
			file4_name.append(output_path);
			file4_name.append("/");
			file4_name.append(lib_name);
			file4_name.append(".group_");
			file4_name.append(group_no);
			file4_name.append(".2.");
			file4_name.append(file_format);
			file4_name.append(".gz");

			copyfile(file1_name,file3_name);
			copyfile(file2_name,file4_name);
			groups_base_num[final_group_st[0]] += groups_base_num[final_group_st[id]];
			groups_base_num[final_group_st[id]] = 0;
		}

	}else
	{
		sprintf(temp_group,"%d",final_group_st[0]);
		group_no = (string)temp_group;

		string file1_name;
		file1_name.reserve(output_path.size() + 50);
		file1_name.append(output_path);
		file1_name.append("/");
		file1_name.append(lib_name);
		file1_name.append(".group_");
		file1_name.append(group_no);
		file1_name.append(".1.");
		file1_name.append(file_format);
		file1_name.append(".gz");
				
		for(vector<size_t>::size_type id = 1;id != final_group_st.size();++id)
		{
			sprintf(temp_group,"%d",final_group_st[id]);
			group_no = (string)temp_group;

			string file3_name;
			file3_name.reserve(output_path.size() + 50);
			file3_name.append(output_path);
			file3_name.append("/");
			file3_name.append(lib_name);
			file3_name.append(".group_");
			file3_name.append(group_no);
			file3_name.append(".1.");
			file3_name.append(file_format);
			file3_name.append(".gz");

			copyfile(file1_name,file3_name);
			groups_base_num[final_group_st[0]] += groups_base_num[final_group_st[id]];
			groups_base_num[final_group_st[id]] = 0;

		}
	}
}

//reads输出到文件,追加方式
void output_reads_to_file(string &lib_name,string group_no,string &file_format,string file_no,
			  string &seq_id,string &seq_base,string &seq_plus,string &seq_qual)
{
	string file_name;
	file_name.reserve(output_path.size() + 50);

	file_name.append(output_path);
	file_name.append("/");
	file_name.append(lib_name);
	if(!file_format.compare("fastq"))
	{	
		file_name.append(".group_");
		file_name.append(group_no);
		file_name.append(".");
		file_name.append(file_no);
		file_name.append(".");
		file_name.append(file_format);
		file_name.append(".gz");
       		ogzstream file(file_name.c_str(),std::ios::app);
	        file << seq_id << endl;
		file << seq_base << endl;
		file << seq_plus << endl;
		file << seq_qual << endl;
	        file.close();
	}else
	{
		file_name.append(".group_");
		file_name.append(group_no);
		file_name.append(".");
		file_name.append(file_no);
		file_name.append(".");
		file_name.append(file_format);
		file_name.append(".gz");
	       	ogzstream file(file_name.c_str(),std::ios::app);
		file << ">" << seq_id << "\n" << seq_base;		//seq_base,包含和换行符
        	file.close();
	}
}

//将reads的uniq kmer指向指定的group
void reads_assign_to_group(vector<Kmer *> *uniq_kmer_address,unsigned int group_id)
{
	for(vector<Kmer *>::size_type idx = 0;idx != uniq_kmer_address->size();++idx)
	{
		bool flag = false;//true说明已经包含该group
		GList * temp_group = uniq_kmer_address->at(idx)->group;

		GList *last_group = NULL;
		while(temp_group)
		{
			if(temp_group->group_id < group_id)
			{
				last_group = temp_group;
				temp_group = temp_group->next;
			}
			else if(temp_group->group_id == group_id)
			{
				flag = true;
				break;	
			}else
			{
				break;
			}
		}
		
		if(flag)continue;

		//在链表头部插入;一种是尚未有链表，另一种是group_id 小于表头
		GList *new_group = new GList(group_id);
		new_group->next = temp_group;
		if(last_group)//在中部和尾部插入
		{
			last_group->next = new_group;
		}else//在头部插入
		{
			uniq_kmer_address->at(idx)->group = new_group;
		}
	}
}

//多个group的情况，刷新tree的Group信息
//第一个group保留，其他的group整合到第一个group上
void flush_group_infor()
{
	//启用多线程
	pthread_mutex_init(&flush_lock,NULL);
	hash_idx = 0;
	pthread_t *tids = new pthread_t[flush_thread_num];
	for(int i = 0;i != flush_thread_num ;++i)
	{
		int ret = pthread_create(&tids[i],NULL,tranvers_hash,NULL);
		//if(DEBUG){cout << "create flush  thread " <<tids[i] << endl;}
		if(ret != 0)
		{
			cerr << "pthread creat error: error code " << ret << endl;
		}
	}
	//等待进程
	for(int i = 0;i != flush_thread_num;++i)
	{
		pthread_join(tids[i],NULL);
		//cout << tids[i] << " exit.flush" << endl;
	}

	//group_infor 置零
	for(vector<size_t>::size_type id = 1; id != final_group_st.size();++id)
	{
		group_infor[final_group_st[id]] = false;
	}
	delete[] tids;
	tids = NULL;
}

//根据uniq_kmer数量、uniq_kmer_address、min_kmer_num、hash完成reads的分配。
void select_groups(vector<Kmer *> *uniq_kmer_address,int min_kmer_num)
{
	int group_kmer_num[MAX_GROUP_NUM];					//统计每个group的uniq kmer数
	memset(group_kmer_num,0,sizeof(group_kmer_num));
	final_group_st.clear();
	final_group_ui.clear();
	for(vector<Kmer *>::size_type idx = 0;idx != uniq_kmer_address->size();++idx)
	{
		GList * temp_group = uniq_kmer_address->at(idx)->group;
		while(temp_group != NULL)
		{
			++group_kmer_num[temp_group->group_id];
			temp_group = temp_group->next;
		}

	}

	unsigned int ui_id = 0;
	for(size_t id = 0;id != group_num_id;++id,++ui_id)
	{
		//if(group_kmer_num[id] >= min_kmer_num){final_group_st.push_back(id);final_group_ui.push_back(ui_id);}
		if(group_kmer_num[id] >= 1){final_group_st.push_back(id);final_group_ui.push_back(ui_id);}
	}
}

bool find_uninitial_group(unsigned int &uninitial_id)
{
	unsigned int gid = 0;
	for(size_t id = 0;id != group_infor.size();++id,++gid)
	{
		if(!group_infor[id])
		{
			uninitial_id = gid;
			return true;
		}
	}
	return false;
}

//判断该lib是否已经存在数组中
bool is_not_exit_lib(string &str,vector<string> &vec)
{
	for(vector<string>::size_type id= 0; id != vec.size(); ++id)
	{
		if(!str.compare(vec[id]))return false;
	}
	return true;
}

//分PE reads,分配成功返回true
bool is_finish_assigned_pe_reads(int reads_type,int uniq_kmer_num,int min_kmer_num,int initialize_kmer_num,
				 vector<Kmer *> *uniq_kmer_address,string &file_format,string &library_name ,bool last_operate,
				string &reads1_id,string &reads1_seq,string &reads1_plus,string &reads1_qual,
				string &reads2_id,string &reads2_seq,string &reads2_plus,string &reads2_qual,long reads_base_num)
{
	//if(uniq_kmer_num <= min_kmer_num && !last_operate)return false;		//uniq_kmer总数小于阈值，直接退出
	if(uniq_kmer_num <= 0 && !last_operate)return false;		//uniq_kmer总数小于阈值，直接退出 //debug

	//这里开始加分配锁
	pthread_mutex_lock(&assign_lock);

	//更新final_group_st
	select_groups(uniq_kmer_address,min_kmer_num);

	vector<size_t>::size_type group_num = final_group_st.size();
	char temp_group[5];
	string group_no;
	if(group_num == 0)					//不属于任何一个group
	{
		if(last_operate)
		{
			unsigned int group_id;
			if(uniq_kmer_num >= initialize_kmer_num && find_uninitial_group(group_id))//大于初始化的值,可以初始化
			{
				sprintf(temp_group,"%d",group_id);
				group_no = (string)temp_group;

				reads_assign_to_group(uniq_kmer_address,group_id);
				output_reads_to_file(library_name,group_no,file_format,"1",reads1_id,reads1_seq,reads1_plus,reads1_qual);
				if(reads_type == 2)
				{
					output_reads_to_file(library_name,group_no,file_format,"2",reads2_id,reads2_seq,reads2_plus,reads2_qual);
				}
				group_infor[group_id] = true;
				groups_base_num[group_id] = reads_base_num;
			}else//不能初始化group，输出到remain
			{
				output_reads_to_file(library_name,"remain",file_format,"1",reads1_id,reads1_seq,reads1_plus,reads1_qual);
				if(reads_type == 2)
				{
					output_reads_to_file(library_name,"remain",file_format,"2",reads2_id,reads2_seq,reads2_plus,reads2_qual);
				}
				groups_base_num[MAX_GROUP_NUM] += reads_base_num;
				

			}
		}else
		{
			pthread_mutex_unlock(&assign_lock);
			return false;	
		}
	}
	else if(group_num == 1)					//属于一个group
	{
		sprintf(temp_group,"%d",final_group_st[0]);
		group_no = (string)temp_group;
		reads_assign_to_group(uniq_kmer_address,final_group_st[0]);
		output_reads_to_file(library_name,group_no,file_format,"1",reads1_id,reads1_seq,reads1_plus,reads1_qual);
		if(reads_type == 2)
		{
			output_reads_to_file(library_name,group_no,file_format,"2",reads2_id,reads2_seq,reads2_plus,reads2_qual);
		}
		groups_base_num[final_group_st[0]] += reads_base_num;
	}
	else							//属于多个group
	{
		sprintf(temp_group,"%d",final_group_st[0]);
		group_no = (string)temp_group;
		reads_assign_to_group(uniq_kmer_address,final_group_st[0]);
		flush_group_infor();			//递归刷新group信息
		merge_file(reads_type,library_name,file_format);	//文件合并
		output_reads_to_file(library_name,group_no,file_format,"1",reads1_id,reads1_seq,reads1_plus,reads1_qual);
		if(reads_type == 2)
		{
			output_reads_to_file(library_name,group_no,file_format,"2",reads2_id,reads2_seq,reads2_plus,reads2_qual);
		}
		groups_base_num[final_group_st[0]] += reads_base_num;

	}
	pthread_mutex_unlock(&assign_lock);
	return true;
}

//分一个PE文件的reads
void* classify_one_lib(void *ar)
{
	finish_reads_num = 0;
	int reads_type = atoi(current_file[0].c_str());			//1为SE，2为PE
	string library_name = current_file[1];				//文库编号。相同文库会合并
	string file_format = current_file[2];				//文件类型，fasta/fastq
	long reads_base_num = 0;

	string reads1_id,reads1_seq,reads1_plus,reads1_qual;
	string reads2_id,reads2_seq,reads2_plus,reads2_qual;
	vector<Kmer *> *uniq_kmer_address = new vector<Kmer *>;
	uniq_kmer_address->reserve(200000);


	pthread_t pid = pthread_self();
	while(true)
	{
		reads1_seq.clear();
		reads2_seq.clear();
		pthread_mutex_lock(&reads_lock);
		if(!reads_file1.eof())
		{
			++finish_reads_num;
			if(finish_reads_num % 10000 == 0)cout << "finish " << finish_reads_num << " reads." << endl;
			if(!file_format.compare("fastq"))
			{
				getline(reads_file1,reads1_id);
				getline(reads_file1,reads1_seq);
				getline(reads_file1,reads1_plus);
				getline(reads_file1,reads1_qual);
				if(reads_type == 2)
				{
					getline(reads_file2,reads2_id);
					getline(reads_file2,reads2_seq);
					getline(reads_file2,reads2_plus);
					getline(reads_file2,reads2_qual);
				}
			}else if(!file_format.compare("fasta"))
			{
				getline(reads_file1,reads1_id);
				getline(reads_file1,reads1_seq,'>');
				if(reads_type == 2)
				{
					getline(reads_file2,reads2_id);
					getline(reads_file2,reads2_seq,'>');
				}
			}else
			{
				cerr << "reads file format must fastq or fasta." <<endl;
				exit(0);
			}

			pthread_mutex_unlock(&reads_lock);
			
			//将序列转化成正确的格式，除ATCGN以外的字符都会删除
			string new_reads1_seq;
			string new_reads2_seq;
			new_reads1_seq.reserve(reads1_seq.size());
			string_clear(reads1_seq,new_reads1_seq);
			reads_base_num = new_reads1_seq.length();
			int reads_kmer_num = reads_base_num - KMER_SIZE + 1;
			if(reads_kmer_num < 0)continue;
			if(reads_type == 2)
			{
				new_reads2_seq.reserve(reads2_seq.size());
				string_clear(reads2_seq,new_reads2_seq);
				reads_kmer_num += reads_kmer_num;//翻倍
				reads_base_num += reads_base_num;
			}
			 
			int min_kmer_num = (int)(reads_kmer_num * overlap_threshold);
			int initialize_kmer_num = (int)(reads_kmer_num * initialize_threshold);
			int uniq_kmer_num = 0;					//uniq_kmer 的个数


			if(reads_type == 2)
			{
				//uniq_kmer_address->clear(); //stat_pe_uniq_kmer中有clear操作
				uniq_kmer_num = stat_pe_uniq_kmer(new_reads1_seq,new_reads2_seq,uniq_kmer_address);//reads2 取反链
			}else 
			{
				uniq_kmer_address->clear(); //stat_pe_uniq_kmer中有clear操作
				uniq_kmer_num = stat_reads_uniq_kmer(new_reads1_seq,uniq_kmer_address,true);
			}
		

			if(is_finish_assigned_pe_reads(reads_type,uniq_kmer_num,min_kmer_num,initialize_kmer_num,uniq_kmer_address,
			   file_format,library_name,false,reads1_id,reads1_seq,reads1_plus,reads1_qual,
			   reads2_id,reads2_seq,reads2_plus,reads2_qual,reads_base_num))continue;

			//------------------------------------互补链----------------------------------------------------
			if(reads_type == 2)
			{
				//uniq_kmer_address->clear(); //stat_pe_uniq_kmer中有clear操作
				uniq_kmer_num = stat_pe_uniq_kmer(new_reads2_seq,new_reads1_seq,uniq_kmer_address);//reads1 取反链
			}else 
			{
				uniq_kmer_address->clear(); //stat_pe_uniq_kmer中有clear操作
				uniq_kmer_num = stat_reads_uniq_kmer(new_reads1_seq,uniq_kmer_address,false);
			
}
			//把uniq kmer数少于阈值的扔到垃圾桶
			if(uniq_kmer_num < min_kmer_num)
			{
				output_reads_to_file(library_name,"garbage",file_format,"1",reads1_id,reads1_seq,reads1_plus,reads1_qual);
				if(reads_type == 2)
				{
					output_reads_to_file(library_name,"garbage",file_format,"2",reads2_id,reads2_seq,reads2_plus,reads2_qual);
				}
				groups_base_num[MAX_GROUP_NUM] += reads_base_num;
				continue;
			}

			is_finish_assigned_pe_reads(reads_type,uniq_kmer_num,min_kmer_num,initialize_kmer_num,uniq_kmer_address,
			   file_format,library_name,true,reads1_id,reads1_seq,reads1_plus,reads1_qual,
			   reads2_id,reads2_seq,reads2_plus,reads2_qual,reads_base_num);


		}else
		{
			pthread_mutex_unlock(&reads_lock);
			break;
		}
	}
	cout << "thread " << pid <<" exit." << endl;
}

//分reads
void classify_reads(vector<vector<string> > &all_reads_infor,int iterate_times)
{
	vector<Kmer *> *uniq_kmer_address = new vector<Kmer *>;		//为兼容大的contigs初始化，使用堆内存
	uniq_kmer_address->reserve(20000);				//uniq_kmer_address 设定为最多的情况
	int reads_type;
	string library_name,file_format,reads1_file_name,reads2_file_name;	

	//用于迭代
	vector<string> added_library;
	added_library.reserve(10);
	int it = 0;
	vector<string>::size_type size_flag = all_reads_infor.size() - 1;

	for(vector<string>::size_type index = 0;index != all_reads_infor.size();++index)
	{
		current_file = all_reads_infor[index];
		reads_type = atoi(current_file[0].c_str());			//1为SE，2为PE
		library_name = current_file[1];				//文库编号。相同文库会合并
		file_format = current_file[2];				//文件类型，fasta/fastq
		reads1_file_name = current_file[3];				//reads1文件名

		cout << "*********************************************************" <<endl;
				
		if(reads_type == 2)							//PE reads
		{
			reads2_file_name = current_file[4];

                        //改变remain的文件名，因为remain的文件名需要作为输出文件名
                        if(it != 0)
                        {
                                string temp_file1 = reads1_file_name;
                                string temp_file2 = reads2_file_name;

                                temp_file1.append(".temp.gz");
                                temp_file2.append(".temp.gz");

                                rename(reads1_file_name.c_str(),temp_file1.c_str());
                                rename(reads2_file_name.c_str(),temp_file2.c_str());

                                reads1_file_name = temp_file1;
                                reads2_file_name = temp_file2;
                        }
			reads_file1.open(reads1_file_name.c_str());
			reads_file2.open(reads2_file_name.c_str());
			if(!reads_file1){cerr << "fail to open input file : " << reads1_file_name << endl;exit(1);}
			if(!reads_file2){cerr << "fail to open input file : " << reads2_file_name << endl;exit(1);}

			//启用多线程
			pthread_mutex_init(&reads_lock,NULL);
			pthread_mutex_init(&assign_lock,NULL);
			pthread_t *tids = new pthread_t[calculate_thread_num];
			for(int i = 0;i != calculate_thread_num ;++i)
			{
				int ret = pthread_create(&tids[i],NULL,classify_one_lib,NULL);
				if(DEBUG)cout << "create thread " <<tids[i] << endl;
				if(ret != 0)
				{
					cerr << "pthread creat error: error code " << ret << endl;
				}
			}

			//等待进程
			for(int i = 0;i != calculate_thread_num;++i)
			{
				pthread_join(tids[i],NULL);
			}
			delete []tids;
			tids=NULL;

			reads_file1.close();
			reads_file1.clear();
			reads_file2.close();
			reads_file2.clear();

                         //删除临时文件
                        if(it != 0)
                        {
                              remove(reads1_file_name.c_str());
                              remove(reads2_file_name.c_str());
                        }


		}else if(reads_type == 1)					//SE reads
		{
			reads2_file_name = "";
			//改变remain的文件名，因为remain的文件名需要作为输出文件名
                        if(it != 0)
                        {
                                string temp_file1 = reads1_file_name;
                                temp_file1.append(".temp.gz");
                                rename(reads1_file_name.c_str(),temp_file1.c_str());
                                reads1_file_name = temp_file1;
                        }

			reads_file1.open(reads1_file_name.c_str());
			if(!reads_file1){cerr << "fail to open input file : " << reads1_file_name << endl;exit(1);}
		
			//启用多线程
			pthread_attr_t attr;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);
			pthread_mutex_init(&reads_lock,NULL);
			pthread_mutex_init(&assign_lock,NULL);
			pthread_t *tids = new pthread_t[calculate_thread_num];
			for(int i = 0;i != calculate_thread_num ;++i)
			{
				int ret = pthread_create(&tids[i],NULL,classify_one_lib,NULL);
				if(DEBUG)cout << "create thread " <<tids[i] << endl;
				if(ret != 0)
				{
					cerr << "pthread creat error: error code " << ret << endl;
				}
			}
			//等待进程
			for(int i = 0;i != calculate_thread_num;++i)
			{
				pthread_join(tids[i],NULL);
			}
			delete []tids;
			tids=NULL;
			pthread_attr_destroy(&attr);

			reads_file1.close();
			reads_file1.clear();

                        //删除临时文件
                        if(it != 0)
                        {
                                remove(reads1_file_name.c_str());
                        }
		}else								//error reads
		{
			cerr << "ERROR: reads type error.must be 1 or 2."<<endl;
			exit(0);
		}
		system("date");
		cout << "finish split " << reads1_file_name << "," << reads2_file_name << endl;

		//当还需要迭代时，且该lib未加入队列，将lib加入队列
		if(iterate_times != 0 && is_not_exit_lib(library_name,added_library))
		{
			string remain_file1;
			string remain_file2;
			remain_file1.reserve(output_path.size() + 50);
			remain_file1.append(output_path);
			remain_file1.append("/");
			remain_file1.append(library_name);
			remain_file1.append(".group_remain.1.");
			remain_file1.append(file_format);
			remain_file1.append(".gz");

			vector<string> reads_file;
			reads_file.reserve(5);
			if(reads_type == 1)
			{
				reads_file.push_back("1");
			}else
			{
				reads_file.push_back("2");
				remain_file2.reserve(output_path.size() + 50);
				remain_file2.append(output_path);
				remain_file2.append("/");
				remain_file2.append(library_name);
				remain_file2.append(".group_remain.2.");
				remain_file2.append(file_format);
				remain_file2.append(".gz");
			}
			reads_file.push_back(library_name);
			reads_file.push_back(file_format);
			reads_file.push_back(remain_file1);
			reads_file.push_back(remain_file2);

			all_reads_infor.push_back(reads_file);
			added_library.push_back(library_name);
		}

		//输入当前group的base信息
		for(size_t gi = 0;gi <= group_num_id;++gi)
		{
			cout << library_name << " group " << gi << ": " << groups_base_num[gi] << "bp." << endl;
		}

		//更新size_flag
		if(iterate_times != 0 && index == size_flag)//上一次迭代完成,
		{
			//memset(groups_base_num,0,sizeof(groups_base_num));//清零
			groups_base_num[group_num_id] = 0;
			added_library.clear();
			size_flag = all_reads_infor.size() - 1;
			--iterate_times;
			++it;
			cout << "#######################################################" <<endl;
			cout << "#                Iteratr "<< it << " times                  #" << endl;
			cout << "#######################################################" <<endl;
		}
	}
	delete uniq_kmer_address;
	uniq_kmer_address = NULL;
}


int main(int argc, char *argv[])
{
	int c;
	string kmer_file_name;
	string reads_file_name;
	initialize_threshold = 0.8;
	overlap_threshold = 0.2;
	int iterate_times = 1;
	while((c=getopt(argc,argv,"k:r:o:I:i:p:c:f:g:")) != -1)
	{
		switch(c)
		{
			case 'k' : kmer_file_name = optarg;break;
			case 'r' : reads_file_name = optarg;break;
			case 'p' : output_path = optarg;break;
			case 'i' : initialize_threshold = atof(optarg);break;
			case 'o' : overlap_threshold = atof(optarg);break;
			case 'I' : iterate_times = atoi(optarg);break;
			case 'c' : calculate_thread_num = atoi(optarg);break;
			case 'f' : flush_thread_num = atoi(optarg);break;
			case 'g' : group_num_no = (unsigned int)atoi(optarg);break;
			default  : cout<<"error:"<<(char)c<<endl;break;
		}
	}
	group_num_id = (size_t)group_num_no;//小变大
	if(group_num_id > MAX_GROUP_NUM){cerr << "group number can not greater than " << MAX_GROUP_NUM << endl; exit(0);}
	//open kmer file
	igzstream kmer_file;
	kmer_file.open(kmer_file_name.c_str());
	if(!kmer_file)
	{
		cerr << "fail to open input file : " << kmer_file_name << endl;
		exit(1);
	}

	ifstream reads_file;
	reads_file.open(reads_file_name.c_str());
	if(!reads_file)
	{
		cerr << "fail to open input file : " << reads_file_name << endl;
		exit(1);
	}

//----------------------------------------	
//读入reads信息
//----------------------------------------	
	vector<vector<string> > all_reads_infor;
	all_reads_infor.reserve(10);
	string reads_infor;
	string deli = "\t";
	vector<string> fragment;
	fragment.reserve(4);				//设定4个预留空间
	while(getline(reads_file,reads_infor))
	{
		fragment.clear();			//清空vector，但没有释放内存
		split_string(reads_infor,fragment,deli,false);
		all_reads_infor.push_back(fragment);
	}
	reads_file.close();
	vector<string>().swap(fragment);		//释放内存
//----------------------------------------
//初始化 kmer_complement 和reads_complement
//----------------------------------------
for(int i = REST_KMER_BI_SIZE;i != 0;--i)kmer_complement[--i] = true;
for(int i = MAX_READS_BI_LENGTH;i != 0;--i)reads_complement[--i] = true;
//print_bitset(kmer_complement);
	
//----------------------------------------	
//构建4叉树
//----------------------------------------	

	cout << "----------------------------" <<endl;
	system("date");
	cout <<"creating 4_branch tree..." <<endl;
	hash = new Kmer*[BUCKET_SIZE];
	memset(hash,0,sizeof(hash));

	//creat 4_branch_tree
	string kmer_seq;
	int position;						//用于insert_to_bitset 递归传递参数
	int index;						//用于insert_to_bitset ，节省栈内存
	bitset<HASH_BI_SIZE> hash_bit;
	while(getline(kmer_file, kmer_seq))
	{
		position = 0;
		Kmer *new_kmer = new Kmer;
		kmer_to_bitset(kmer_seq,new_kmer->k,hash_bit);
		//cout << hash_bit << "\t" << new_kmer->k <<endl;
		unsigned long hash_index = hash_bit.to_ulong();
		if(hash[hash_index])				//该节点已经初始化
		{
			Kmer * kmer_flag = hash[hash_index];
			insert_to(new_kmer->k,new_kmer,kmer_flag,position,index);
		}else						//该节点尚未初始化
		{
			Kmer *head_kmer = new Kmer;
			hash[hash_index] = head_kmer;
			Kmer *kmer_flag = head_kmer;
			insert_to(new_kmer->k,new_kmer,kmer_flag,position,index);
		}
	}
	kmer_file.close();
	system("date");
	cout <<"finish  4_branch tree." <<endl;
	cout << "----------------------------" <<endl;
	//tranvers_hash(hash);
	//finish 4_branch_tree

//----------------------------------------	
//初始化group，选取第一个lib的前若干条reads
//----------------------------------------	
	cout << "----------------------------" <<endl;
	system("date");
	cout << "begin to in initialize group infor..." <<endl;
	//cout << group_infor << endl;
	initialize_group(all_reads_infor[0]);
	//cout << group_infor << endl;
	//tranvers_hash(hash);
	system("date");
	cout << "finish initialize group infor." <<endl;
	cout << "----------------------------" <<endl;

//----------------------------------------	
//开始分reads
//----------------------------------------	

	cout << "----------------------------" <<endl;
	system("date");
	cout << "begin to split reads ..." <<endl;
	memset(groups_base_num,0,sizeof(groups_base_num));
	classify_reads(all_reads_infor,iterate_times);
	pthread_mutex_destroy(&reads_lock);
	pthread_mutex_destroy(&assign_lock);
	system("date");
	cout << "finish split all reads." <<endl;
	cout << "----------------------------" <<endl;

}
