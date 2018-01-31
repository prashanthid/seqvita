#include <iostream>
#include <omp.h>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <set>
#include <sys/time.h>
#include "fisher.hpp"
using namespace std;

//default Parameters

int min_coverage = 8;
int min_reads2 = 2;
int min_avg_qual = 15;
float min_var_freq = 0.01;
float min_freq_for_hom = 0.75;
float p_value = 0.99;
int strand_filter =1;

struct data
{
  int forward_ref;
  int reverse_ref;
  int forward_var;
  int reverse_var;
  int quality_ref;
  int quality_var;
  int quality_depth;
  int var_flag;
  string variant_type = "";
  string variant = "";
  bool is_homo;
  bool is_hetero;

  void setzero()
  {
    forward_ref = 0;
    reverse_ref = 0;
    quality_ref = 0;
    var_flag = 0;
    quality_depth = 0;
    quality_var = 0;
    forward_var = 0;
    reverse_var = 0;
  }
};

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

string itoa(int i)
{
  stringstream ss;
  ss << i;
  string a = ss.str();
  return a;
}

string ftoa(float f)
{
  stringstream ss;
  ss << f;
  string a = ss.str();
  return a;
}

std::vector<std::string> Get_Pvalue(int a, int b)
{
  FILE *fpipe;
  char c[100];
  string ref = itoa(a);
  string var = itoa(b);
  string cmd = "Rscript p-value.R " + ref + " " + var;
  if (0 == (fpipe = (FILE*)popen(cmd.c_str(), "r")))
  {
      perror("popen() failed.");
      exit(1);
  }
  fgets(c,sizeof(c),fpipe);
  string str = string(c);
  str = str.substr(0,str.size()-2);
  vector<string> parts;
  parts = split(str,' ',parts);
  parts[1] = parts[1].substr(1,parts[1].size());
  //std::cout << parts[1] << std::endl;
  double pval;
  int flag = 0;
  try {pval = std::stod(parts[1]);}
  catch(std::out_of_range){  for(int i=0;i<parts[1].size();i++)
    {
      if(flag == 0)
      {
        if(i==6)
        {
          if(int(parts[1][i]) > int('4'))
          {
            parts[1][i-1]++;
          }
          parts[1].erase(i,1);
          i--;
          flag = 1;
        }
      }
      else if(i>5 && parts[1][i]!= 'e')
      {
        parts[1].erase(i,1);
        i--;
      }
      else if(parts[1][i] == 'e')
      {
        if(parts[1][i+2] == '0')
        {
          parts[1].replace(i,i,"E" + parts[1].substr(i+1,parts[1].size()));
          parts[1].erase(i+2,1);
          break;
        }
        else
        {
          parts[1].replace(i,i,"E" + parts[1].substr(i+1,parts[1].size()));
          break;
        }
      }
    }
    pclose(fpipe);
    vector <string> p_value;
    p_value.push_back(parts[1]);
    p_value.push_back(parts[2]);
    return p_value;
  }
  stringstream ss;
  ss << std::scientific;
  ss << pval << std::endl;
  string temp =ss.str();
  //std::cout << temp + "\n" + parts[1] << std::endl;
  parts[1] = temp.substr(0,temp.size()-1);
  for(int i=0;i<parts[1].size();i++)
  {
    if(flag == 0)
    {
      if(i==6)
      {
        if(int(parts[1][i]) > int('4'))
        {
          parts[1][i-1]++;
        }
        parts[1].erase(i,1);
        i--;
        flag = 1;
      }
    }
    else if(i>5 && parts[1][i]!= 'e')
    {
      parts[1].erase(i,1);
      i--;
    }
    else if(parts[1][i] == 'e')
    {
      if(parts[1][i+2] == '0')
      {
        parts[1].replace(i,i,"E" + parts[1].substr(i+1,parts[1].size()));
        parts[1].erase(i+2,1);
        break;
      }
      else
      {
        parts[1].replace(i,i,"E" + parts[1].substr(i+1,parts[1].size()));
        break;
      }
    }
  }
  pclose(fpipe);
  vector <string> p_value;
  /*if(pval == 0){
    p_value.push_back("0");
  }
  else{
    p_value.push_back(parts[1]);
  }*/
  p_value.push_back(parts[1]);
  p_value.push_back(parts[2]);
  return p_value;
}

void vcf_Header(int minAvgQual,int no_of_samples)
{
  string vcf_header = "##fileformat=VCFv4.1\n";
  vcf_header += "##source=VarScan2\n";
  vcf_header += "##INFO=<ID=ADP,Number=1,Type=Integer,Description=\"Average per-sample depth of bases with Phred score >= " + itoa(minAvgQual) + "\">\n";
  vcf_header += "##INFO=<ID=WT,Number=1,Type=Integer,Description=\"Number of samples called reference (wild-type)\">\n";
  vcf_header += "##INFO=<ID=HET,Number=1,Type=Integer,Description=\"Number of samples called heterozygous-variant\">\n";
  vcf_header += "##INFO=<ID=HOM,Number=1,Type=Integer,Description=\"Number of samples called homozygous-variant\">\n";
  vcf_header += "##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Number of samples not called\">\n";
  vcf_header += "##FILTER=<ID=str10,Description=\"Less than 10% or more than 90% of variant supporting reads on one strand\">\n";
  vcf_header += "##FILTER=<ID=indelError,Description=\"Likely artifact due to indel reads at this position\">\n";
  vcf_header += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
  vcf_header += "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
  vcf_header += "##FORMAT=<ID=SDP,Number=1,Type=Integer,Description=\"Raw Read Depth as reported by SAMtools\">\n";
  vcf_header += "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Quality Read Depth of bases with Phred score >= " + itoa(minAvgQual) + "\">\n";
  vcf_header += "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">\n";
  vcf_header += "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">\n";
  vcf_header += "##FORMAT=<ID=FREQ,Number=1,Type=String,Description=\"Variant allele frequency\">\n";
  vcf_header += "##FORMAT=<ID=PVAL,Number=1,Type=String,Description=\"P-value from Fisher's Exact Test\">\n";
  vcf_header += "##FORMAT=<ID=RBQ,Number=1,Type=Integer,Description=\"Average quality of reference-supporting bases (qual1)\">\n";
  vcf_header += "##FORMAT=<ID=ABQ,Number=1,Type=Integer,Description=\"Average quality of variant-supporting bases (qual2)\">\n";
  vcf_header += "##FORMAT=<ID=RDF,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on forward strand (reads1plus)\">\n";
  vcf_header += "##FORMAT=<ID=RDR,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases on reverse strand (reads1minus)\">\n";
  vcf_header += "##FORMAT=<ID=ADF,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on forward strand (reads2plus)\">\n";
  vcf_header += "##FORMAT=<ID=ADR,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases on reverse strand (reads2minus)\">\n";
  vcf_header +=  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  for(int colCounter = 0; colCounter < no_of_samples; colCounter++)
  {
    vcf_header += "\tSample" + itoa(colCounter+1);
  }
  //std::cout << vcf_header << std::endl;
  printf("%s\n",vcf_header.c_str());
}

void pileuptovcf(std::vector<string> &pileup, std::vector<data>  &individual)//convert the variant obtained into vcf format
{
  //std::cout << "Entered" << '\n';
  //std::cout << individual[0].quality_depth << '\n';
  int i;
  int no_of_samples = individual.size();
  char Reference = pileup[2][0];
  if(pileup[2][0] <= 'z' && pileup[2][0] >= 'a')
  {
    Reference = pileup[2][0] - 'a' + 'A';
  }
  string vcf_line = pileup[0] + '\t' + pileup[1] + '\t' + '.' + '\t';
  int position = 0;
  int deletion_size = 0;
  string ref(1,Reference);
  for (i=0;i<no_of_samples;i++){
    if(individual[i].var_flag && individual[i].variant_type == "DEL"){
        deletion_size = individual[i].variant.length();
        if(deletion_size > ref.length()-1){
          ref = Reference;
          ref += individual[i].variant;
      }
    }
  }
  vcf_line += ref + '\t';
  set <string> variant_list;
  map <string,int> variant_order;
  int variant_number = 1;
  for(i=0;i<no_of_samples;i++){
    if(individual[i].var_flag){
      if(individual[i].variant_type == "DEL"){
        string del = ref;
        del.erase(ref.length()-individual[i].variant.length());
        individual[i].variant.clear();
        individual[i].variant += del;
      }
      if(individual[i].variant_type == "INS"){
        string ins(1,ref[0]);
        ins += individual[i].variant + ref.substr(1);
        individual[i].variant.clear();
        individual[i].variant += ins;
      }
      if(individual[i].variant_type == "SNP"){
        string snp = ref;
        snp[0] = individual[i].variant[0];
        individual[i].variant.clear();
        individual[i].variant += snp;
      }
      if(variant_list.find(individual[i].variant) == variant_list.end()){
        variant_list.insert(individual[i].variant);
        variant_order[individual[i].variant] = variant_number++;
        vcf_line += individual[i].variant + ",";
      }
    }
  }
  vcf_line = vcf_line.substr(0,vcf_line.size()-1);
  vcf_line += "\t.\tPASS\t";
  int total_depth = 0;
  int wild_type = 0;
  int homozygous = 0;
  int heterozygous = 0;
  string Genotype = "";
  for(i = 0;i < no_of_samples;i++)
  {
    int ref_depth = individual[i].forward_ref + individual[i].reverse_ref;
    int var_depth = individual[i].forward_var + individual[i].reverse_var;
    //std::cout << ref_depth << '\n';

    Genotype += "\t";
    if(individual[i].var_flag)
    {
      if(individual[i].is_homo)
      {
        homozygous++;
        int geno = variant_order[individual[i].variant];
        Genotype += to_string(geno) + "/" + to_string(geno);
      }
      else if(individual[i].is_hetero)
      {
        heterozygous++;
        int geno = variant_order[individual[i].variant];
        Genotype +=  "0/" + to_string(geno);
      }
      //vector <string> GQ = Get_Pvalue(ref_depth,var_depth);
      int d = (ref_depth+var_depth)*0.001;
      double var_pval = fisher_test(ref_depth,var_depth,individual[i].quality_depth,d);
      stringstream ss;
      ss << std::scientific;
      ss << var_pval << std::endl;
      string pvalstr =ss.str();
      /*if(GQ[1][0] == 'N')
        GQ[1] = ".";
      else{
        int Genotype_Quality = stoi(GQ[1]);
        int GQ_threshold = -10 * log10(p_value);
        //std::cout << Genotype_Quality << '\n';
        if(Genotype_Quality < GQ_threshold)
        {
          individual[i].var_flag = 0;
          i--;
          continue;
        }
        if(Genotype_Quality > 255){
          GQ[1] = "255";
        }
      }*/
      if(var_pval > p_value){
        individual[i].var_flag = 0;
        i--;
        continue;
      }
      int GQ = -10*log10(var_pval);

      Genotype += ":" + itoa(GQ) + pileup[3] + ":" + itoa(individual[i].quality_depth) + ":" + itoa(ref_depth) + ":";
      //int depth = individual[i].forward_var+individual[i].reverse_var;
      Genotype += itoa(var_depth) + ":" ;
      float temp = roundf((float(var_depth)/float(individual[i].quality_depth))*10000);
      Genotype += ftoa(float(temp)/100) + "%:" + pvalstr + ":";
      if(ref_depth)
      {
        Genotype += itoa(individual[i].quality_ref/ref_depth) + ":" + itoa(individual[i].quality_var/var_depth);
      }
      else
      {
        Genotype += string("0") + ":" + itoa(individual[i].quality_var/var_depth);
      }
      Genotype += ":" + itoa(individual[i].forward_ref) + ":" + itoa(individual[i].reverse_ref) + ":" + itoa(individual[i].forward_var) + ":" + itoa(individual[i].reverse_var);
    }
    else
    {
      wild_type++;
      Genotype += "./.:.:.:.:.:.:.:.:.:.:.:.:.:.";
    }
    total_depth += individual[i].quality_depth;
  }
  vcf_line += "ADP=" + itoa(total_depth/no_of_samples) + ";WT=" + itoa(wild_type) + ";HET=" + itoa(heterozygous) + ";HOM=" + itoa(homozygous) + ";NC=" + itoa(no_of_samples - heterozygous - homozygous - wild_type);
  vcf_line += "\tGT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR" + Genotype;
  //std::cout << vcf_line << std::endl;
  printf("%s\n",vcf_line.c_str() );
}

void detectvariant(data &individual,string bases,string qualities,int rawdepth)//detect variant from pileup format
{
  std::map<string, data> variants;
  int prevBase_quality,nextBase_quality;
  int quality_offset=0,quality_depth=0,quality_ref=0;
  int reverse_ref=0,forward_ref=0;
  string nextBase;
  string prevBase;
  //std::cout << "Started" << "\n";
  for(int read_no = 0;read_no < bases.length();read_no++)
  {
    //std::cout << "Entered" << '\n';
    prevBase = "";
    if(read_no > 1 && read_no < bases.length() -1){
      prevBase = bases[read_no-1];
    }

    nextBase = "";
    if(read_no < bases.length()-1){
      nextBase = bases[read_no+1];
    }
    //std::cout << nextBase << '\n';
    if (bases[read_no] == ',' && !(nextBase == "+" || nextBase == "-"))
    {
      char score = qualities[quality_offset++];
      int qual = score-33;
      //std::cout << qual << '\n';
      if(qual >= min_avg_qual)
      {
        quality_depth++;
        reverse_ref += 1;
        quality_ref += qual;
      }
    }
    else if (bases[read_no] == '.' && !(nextBase == "+" || nextBase == "-"))
    {
      char score = qualities[quality_offset++];
      int qual = score-33;
      //std::cout << qual << '\n';
      if(qual >= min_avg_qual)
      {
        quality_depth++;
        forward_ref += 1;
        quality_ref += qual;
      }
    }
    else if(bases[read_no] == 'A' || bases[read_no] == 'T' || bases[read_no] == 'G' || bases[read_no] == 'C')
    {
      char score = qualities[quality_offset++];
      int qual = score-33;
      //std::cout << qual << '\n';
      if(qual >= min_avg_qual)
      {
        string variant(1,bases[read_no]);
        quality_depth++;
        if(variants.find(variant) == variants.end()){
          variants[variant].setzero();
          variants[variant].forward_var = 1;
          variants[variant].quality_var = qual;
          variants[variant].variant_type = "SNP";
          variants[variant].variant = bases[read_no];
        }
        else{
          variants[variant].forward_var += 1;
          variants[variant].quality_var += qual;
        }
      }
    }
    else if(bases[read_no] == 'a' || bases[read_no] == 't' || bases[read_no] == 'g' || bases[read_no] == 'c')
    {
      char score = qualities[quality_offset++];
      int qual = score-33;
      //std::cout << qual << '\n';
      if(qual >= min_avg_qual)
      {
        quality_depth++;
        string variant(1,toupper(bases[read_no]));
        if(variants.find(variant) == variants.end()){
          variants[variant].setzero();
          variants[variant].reverse_var = 1;
          variants[variant].quality_var = qual;
          variants[variant].variant_type = "SNP";
          variants[variant].variant = variant;
        }
        else{
          variants[variant].reverse_var += 1;
          variants[variant].quality_var += qual;
        }
      }
    }
    else if(bases[read_no] == '^')
    {
      read_no++;
    }
    else if(bases[read_no] == '$')
    {
    }
    else if (bases[read_no] == '+')
    {
      char score = qualities[quality_offset++];
      int qual = score-33;
      read_no++;
      string c(1,bases[read_no++]);
      while(bases[read_no] >= '0' && bases[read_no] <= '9')
      {
        c = c + bases[read_no++];
      }
      int skip = atoi(c.c_str());
      string indel(1,bases[read_no]);
      for(int k=1;k<skip;k++){
        indel += bases[read_no+k];
      }
      read_no += skip-1;
      //std::cout << qual << '\n';
      if(qual >= min_avg_qual){
        quality_depth++;
        if(isupper(indel[0])){
          //std::cout << indel << '\n';
          if(variants.find("+"+indel) == variants.end()){
            variants["+"+indel].setzero();
            variants["+"+indel].forward_var = 1;
            variants["+"+indel].quality_var = qual;
            variants["+"+indel].variant_type = "INS";
            variants["+"+indel].variant = indel;
          }
          else{
            variants["+"+indel].forward_var += 1;
            variants["+"+indel].quality_var += qual;
          }
        }
        else{
          transform(indel.begin(), indel.end(), indel.begin(), ::toupper);
          //std::cout << indel << '\n';
          if(variants.find("+"+indel) == variants.end()){
            variants["+"+indel].setzero();
            variants["+"+indel].reverse_var = 1;
            variants["+"+indel].quality_var = qual;
            variants["+"+indel].variant_type = "INS";
            variants["+"+indel].variant = indel;
          }
          else{
            variants["+"+indel].reverse_var += 1;
            variants["+"+indel].quality_var += qual;
          }
        }
      }
    }
    else if (bases[read_no] == '-')
    {
      char score = qualities[quality_offset++];
      int qual = score-33;
      //std::cout << qual << '\n';

      read_no++;
      string c(1,bases[read_no++]);
      while(bases[read_no] >= '0' && bases[read_no] <= '9')
      {
        c = c + bases[read_no++];
      }
      int skip = atoi(c.c_str());
      string indel(1,bases[read_no]);
      for(int k=1;k<skip;k++){
        indel += bases[read_no+k];
      }
      read_no += skip-1;
      if(qual >= min_avg_qual){
        quality_depth++;
        if(isupper(indel[0])){
          if(variants.find("-"+indel) == variants.end()){
            variants["-"+indel].forward_var = 1;
            variants["-"+indel].quality_var = qual;
            variants["-"+indel].variant_type = "DEL";
            variants["-"+indel].variant = indel;
          }
          else{
            variants["-"+indel].forward_var += 1;
            variants["-"+indel].quality_var += qual;
          }
        }
        else{
          transform(indel.begin(), indel.end(), indel.begin(), ::toupper);
          if(variants.find("-"+indel) == variants.end()){
            variants["-"+indel].reverse_var = 1;
            variants["-"+indel].quality_var = qual;
            variants["-"+indel].variant_type = "DEL";
            variants["-"+indel].variant = indel;
          }
          else{
            variants["-"+indel].reverse_var += 1;
            variants["-"+indel].quality_var += qual;
          }
        }
      }
    }
  }
  //std::cout << variants.size() << '\n';
  if(variants.size()>0){
    std::map<string, data>::iterator it=variants.begin();
    //std::cout << it->second.variant << '\n';
    individual = it->second;
    for (;it!=variants.end(); ++it){
      //std::cout << it->first << '\n';
      //std::cout << it->second.forward_var+it->second.reverse_var << '\n';
      if((it->second.forward_var+it->second.reverse_var) > (individual.forward_var+individual.reverse_var)){
        individual = it->second;
      }
    }
  }
  //std::cout << "Ended" << '\n';
  individual.forward_ref = forward_ref;
  individual.reverse_ref = reverse_ref;
  individual.quality_ref = quality_ref;
  individual.quality_depth = quality_depth;
  //std::cout << individual.forward_var+individual.reverse_var << '\n';
}

void findsnp(char ** argv)
{
  std::ifstream file(argv[2]);
  std::string str;
  int number_of_SNPs=0;
  bool vcf_header_flag = 0;
  int sstop = 0;
  //omp_set_num_threads(4);
  #pragma omp parallel private(str)
  {
    while(!sstop)
    {
      #pragma omp critical
      {
        if(!(std::getline(file,str))){
          sstop = 1;
          //break;
        }
      }
      if(sstop){
        #pragma omp flush(sstop)
      }
      if(sstop){
        break;
      }

      vector<string> pileup;
      char delimiter = '	';
      pileup = split(str, delimiter,pileup);
      int no_of_samples = (pileup.size() - 3)/3;
      if(vcf_header_flag == 0)
      {
        vcf_header_flag = 1;
        vcf_Header(min_avg_qual,no_of_samples);
      }
      int i;
      vector <data> individual(no_of_samples);
      bool var_flag = 0;
      //std::cout << pileup[1] << '\n';
      for(i = 0;i < no_of_samples;i++)
      {
        int rawdepth = atoi(pileup[3+3*i].c_str());
        string bases = pileup[4+3*i];
        string qualities = pileup[5+3*i];
        individual[i].setzero();
        //std::cout << "Started" << '\n';
        detectvariant(individual[i],bases,qualities,rawdepth);
        //std::cout << "Ended" << '\n';
        if(individual[i].quality_depth < min_coverage)
        {
          continue;
        }
        else
        {
          int Var_allele_count = 0;
          Var_allele_count = individual[i].forward_var + individual[i].reverse_var;
          if(Var_allele_count < min_reads2)
          {
            continue;
          }
          else
          {
            float strandedness = float(individual[i].forward_var)/float(Var_allele_count);
            float temp = float(Var_allele_count)/float(individual[i].quality_depth);
            if(strandedness < 0.1 || strandedness > 0.9 )
            {
              continue;
            }
            if(temp < min_var_freq)
            {
              continue;
            }
            else if(temp < min_freq_for_hom)
            {
              individual[i].is_hetero = 1;
            }
            else
            {
              individual[i].is_homo = 1;
            }
          }
        }
        individual[i].var_flag = 1;
        if(individual[i].variant_type == "SNP"){
          var_flag = 1;
          number_of_SNPs++;
        }
      }
      if(var_flag)
      {
        pileuptovcf(pileup,individual);
      }

    }
  }
}

void findindel(char ** argv)
{
  std::ifstream file(argv[2]);
  std::string str;
  int number_of_SNPs=0;
  bool vcf_header_flag = 0;
  while(std::getline(file,str))
  {
    vector<string> pileup;
    char delimiter = '	';
    pileup = split(str, delimiter,pileup);
    int no_of_samples = (pileup.size() - 3)/3;
    if(vcf_header_flag == 0)
    {
      vcf_header_flag = 1;
      vcf_Header(min_avg_qual,no_of_samples);
    }
    int i;
    vector <data> individual(no_of_samples);
    bool var_flag = 0;
    for(i = 0;i < no_of_samples;i++)
    {
      int rawdepth = atoi(pileup[3+3*i].c_str());
      string bases = pileup[4+3*i];
      string qualities = pileup[5+3*i];
      individual[i].setzero();
      detectvariant(individual[i],bases,qualities,rawdepth);
      //std::cout << individual[i].quality_depth << '\n';
      if(individual[i].quality_depth < min_coverage)
      {
        continue;
      }
      else
      {
        int Var_allele_count = 0;
        Var_allele_count = individual[i].forward_var + individual[i].reverse_var;
        //std::cout << "Var allele :"<<Var_allele_count << '\n';
        if(Var_allele_count < min_reads2)
        {
          continue;
        }
        else
        {
          float strandedness = float(individual[i].forward_var)/float(Var_allele_count);
          float temp = float(Var_allele_count)/float(individual[i].quality_depth);
          //std::cout << strandedness << '\n';
          if(strandedness < 0.1 || strandedness > 0.9 )
          {
            continue;
          }
          //std::cout << temp << '\n';
          if(temp < min_var_freq)
          {
            continue;
          }
          else if(temp < min_freq_for_hom)
          {
            individual[i].is_hetero = 1;
          }
          else
          {
            individual[i].is_homo = 1;
          }
        }
      }
      individual[i].var_flag = 1;
      if(individual[i].variant_type == "INS" ||individual[i].variant_type == "DEL"){
        var_flag = 1;
        number_of_SNPs++;
      }
    }
    if(var_flag)
    {
      pileuptovcf(pileup,individual);
    }
  }
}

void findgermline(char ** argv)
{
  std::ifstream file(argv[2]);
  std::string str;
  int number_of_SNPs=0;
  bool vcf_header_flag = 0;
  while(std::getline(file,str))
  {
    vector<string> pileup;
    char delimiter = '	';
    pileup = split(str, delimiter,pileup);
    //std::cout << pileup.size() << '\n';
    int no_of_samples = (pileup.size() - 3)/3;
    if(vcf_header_flag == 0)
    {
      vcf_header_flag = 1;
      vcf_Header(min_avg_qual,no_of_samples);
    }
    int i;
    vector <data> individual(no_of_samples);
    bool var_flag = 0;
    for(i = 0;i < no_of_samples;i++)
    {
      int rawdepth = atoi(pileup[3+3*i].c_str());
      string bases = pileup[4+3*i];
      string qualities = pileup[5+3*i];
      individual[i].setzero();
      detectvariant(individual[i],bases,qualities,rawdepth);
      //std::cout << "Quality Depth :" << individual[i].quality_depth << '\n';
      if(individual[i].quality_depth < min_coverage)
      {
        continue;
      }
      else
      {
        int Var_allele_count = 0;
        Var_allele_count = individual[i].forward_var + individual[i].reverse_var;
        //std::cout << "Var allele :"<< individual[i].forward_var+individual[i].reverse_var << '\n';
        if(Var_allele_count < min_reads2)
        {
          continue;
        }
        else
        {
          float strandedness = float(individual[i].forward_var)/float(Var_allele_count);
          float temp = float(Var_allele_count)/float(individual[i].quality_depth);
          if(strandedness < 0.1 || strandedness > 0.9 )
          {
            continue;
          }
          if(temp < min_var_freq)
          {
            continue;
          }
          else if(temp < min_freq_for_hom)
          {
            individual[i].is_hetero = 1;
          }
          else
          {
            individual[i].is_homo = 1;
          }
        }
      }
      individual[i].var_flag = 1;
      var_flag = 1;
      number_of_SNPs++;
      //std::cout << pileup[1] << '\n';
    }
    //std::cout << var_flag << '\n';
    if(var_flag)
    {
      pileuptovcf(pileup,individual);
    }
  }
}

int main (int argc, char *argv[])
{
  bool mpileuptosnp = false;
  bool mpileuptoindel = false;
  bool mpileuptogermline = false;
  if(strcmp(argv[1],"SNP") == 0)
  {
    mpileuptosnp = true;
  }
  else if(strcmp(argv[1],"INDEL") == 0)
  {
    mpileuptoindel = true;
  }
  else if(strcmp(argv[1],"BOTH") == 0)
  {
    mpileuptogermline = true;
  }
  for(int i=3;i< argc;)
  {
    if (strcmp(argv[i],"--min-coverage") == 0)
    {
      min_coverage = atoi(argv[i+1]);
      i += 2;
    }
    else if (strcmp(argv[i],"--min-reads2") == 0)
    {
      min_reads2 = atoi(argv[i+1]);
      i += 2;
    }
    else if (strcmp(argv[i],"--min-avg-qual") == 0)
    {
      min_avg_qual = atoi(argv[i+1]);
      i += 2;
    }
    else if (strcmp(argv[i],"--min-var-freq") == 0)
    {
      min_var_freq = atof(argv[i+1]);
      i += 2;
    }
    else if (strcmp(argv[i],"--min-freq-for-hom") == 0)
    {
      min_freq_for_hom = atof(argv[i+1]);
      i += 2;
    }
    else if (strcmp(argv[i],"--p-value") == 0)
    {
      p_value = atof(argv[i+1]);
      i += 2;
    }
    else if (strcmp(argv[i],"--strand-filter") == 0)
    {
      strand_filter = atoi(argv[i+1]);
      i += 2;
    }
    else
    {
      cout << argv[i] << "  Command not found" << endl;
      exit(0);
    }
  }
  int number_of_SNPs=0;
  //vcf_Header(min_avg_qual);
  bool vcf_header_flag = 0;
  if(mpileuptosnp){
    findsnp(argv);
  }
  else if(mpileuptoindel){
    findindel(argv);
  }
  else if (mpileuptogermline){
    findgermline(argv);
  }
  return 0;
}
