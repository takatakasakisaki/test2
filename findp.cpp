#include <iostream>
#include <execution>
#include <algorithm>
#include <execution>
#include <string>
#include <cmath>
#include <vector>
#include <stdio.h>
#include <string.h>
#include <chrono>
#include <omp.h>
using namespace std;
void find_se(char *line, int *startpos, int *endpos)
{
	int s = -1;
	int e = -1;
	char *p0 = line;
	//find start
	int poss=0;
//printf("line=%s,", line);
	do{
		int c= *p0++;
		if(c == '\0' || c == '\n' || c == '\r'){
			break;
		}
//printf("[%d,%c]", poss, c);
		if(c == '1'){
			if(s < 0){
				s =  poss; //start
			}
			else{
				 e = -1;
			}
		}
		else if(c == '0'){
			if(s >= 0 && e < 0){
				e = poss-1;
			}
		}
		poss ++;
	}while(1);
	if(e < 0){
		e = poss-1;
	}
	//printf(",s=%d,e=%d\n", s, e);
	*startpos = s;
	*endpos = e;

}
int getstartstop(char *fname, std::vector<std::pair<int,int>> &poslist)
{
	int len=0;
	FILE *fp = fopen(fname, "rt");
	if(!fp){
		perror("open");
		return -1;
	}
	do{
		char buf[1024];
		char *p = fgets(buf, sizeof(buf), fp);
		if(!p){
			break;
		}
		int startpos, endpos;
		find_se(buf, &startpos, &endpos);
//printf("%d,%d\n", startpos, endpos);
		std::pair<int,int> startend(startpos,endpos);
		poslist.push_back(startend);
		len++;
	}while(1);
	fclose(fp);
	return len;
}
int main(int argc, char *argv[])
{
	auto ts0 = chrono::system_clock::now();
	auto ts1 = chrono::system_clock::now();
	auto tdiff = chrono::duration_cast<std::chrono::microseconds>(ts1 - ts0);
	cout << tdiff.count();
	char *fname = (char*)"route.txt";
	if(argc < 2){
		printf("fname=%s\n", fname);
	}
	else{
		fname = argv[1];
	}
	printf("fname=%s\n", fname);
	std::vector<std::pair<int,int>> poslist;
	int rv = getstartstop(fname, poslist);
printf("rv=%d\n", rv);
	if(rv > 0){
		for(auto pp : poslist){
			printf("%d,%d\n", pp.first, pp.second);
		}
	}
#if 0
	FILE *fp = fopen(fname, "rt");
	if(!fp){
		perror("open");
		return 1;
	}
	do{
		char buf[1024];
		char *p = fgets(buf, sizeof(buf), fp);
		if(!p){
			break;
		}
		int startpos, endpos;
		find_se(buf, &startpos, &endpos);
	}while(1);
	fclose(fp);
#endif
	return 0;

}

