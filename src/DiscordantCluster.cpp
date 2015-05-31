#include "SnowTools/DiscordantCluster.h"

#include <cassert>

namespace SnowTools {
  
  void DiscordantCluster::addRead(std::string name) {
    std::unordered_map<std::string, bool>::iterator ff = qnames.find(name);
    if (ff == qnames.end())
      qnames.insert(std::pair<std::string, bool>(name, true));
    return;
  }
  
  DiscordantCluster::DiscordantCluster(BamReadVector &this_reads, BamReadVector &all_reads) {
    
    if (this_reads.size() == 0)
      return;
    if (all_reads.size() == 0)
      return;
    
    // check the orientations, fill the reads
    bool rev = this_reads[0].ReverseFlag();
    bool mrev = this_reads[0].MateReverseFlag();
    
    //bool rev = r_is_rstrand(this_reads[0]); //this_reads[0]->IsReverseStrand();
    //bool mrev = r_is_mrstrand(this_reads[0]); //this_reads[0]->IsMateReverseStrand();
    
    //id = r_qname(this_reads[0]); //this_reads[0]->Name;
    id = this_reads[0].Qname();
    
    //cout << "rev " << rev << " mrev " << mrev << endl;
    for (auto& i : this_reads) 
      {
	
	//cout << "flag " << r_flag(i) << " frev " << r_is_rstrand(i) << " mrev " << r_is_mrstrand(i) << endl;
	//assert(rev == r_is_rstrand(i) && mrev == r_is_mrstrand(i)); //i->IsReverseStrand() && mrev == i->IsMateReverseStrand());
	std::string tmp = i.GetZTag("SR");
	//r_get_Z_tag(i, "SR", tmp);
	
	assert(id.length());
	i.AddZTag("DC",id);
	//r_add_Z_tag(i, "DC", id);
	
	assert(tmp.length());
	reads[tmp] = i;
	std::string qn = i.Qname();
	qnames[qn] = true;
	if (qn < id)
	  id = qn;
	if (tmp.at(0) == 't')
	  ++tcount;
	else
	  ++ncount;
      }
    
    addMateReads(all_reads);
    
    assert(reads.size() == mates.size());
    assert(reads.size() > 0);
    
    // set the mapq
    reads_mapq = getMeanMapq(false);
    mates_mapq = getMeanMapq(true);
    
    // set the regions
    reg1 = SnowTools::GenomicRegion(-1,500000000,-1); // read region
    for (auto& i : reads) 
      {
	reg1.strand = i.second.ReverseFlag() ? '-' : '+'; //r_strand(i.second) == '+'; //(!i.second->IsReverseStrand()) ? '+' : '-';
	reg1.chr = i.second.ChrID(); //r_id(i.second); //i.second->RefID;
	if (i.second.Position() < reg1.pos1)
	  reg1.pos1 = i.second.Position(); //r_pos(i.second); //i.second->Position;
	int endpos = i.second.PositionEnd(); //r_endpos(i.second);
	if (endpos > reg1.pos2)
	  reg1.pos2 = endpos;
      }
    
    reg2 = SnowTools::GenomicRegion(-1,500000000,-1); // mate region
    for (auto& i : mates) 
      {
	reg2.strand = i.second.ReverseFlag() ? '-' : '+'; //r_strand(i.second) == '+'; //(!i.second->IsReverseStrand()) ? '+' : '-';
	reg2.chr = i.second.ChrID(); //i.second->RefID;
	if (i.second.Position() < reg2.pos1)
	  reg2.pos1 = i.second.Position();
	int endpos = i.second.PositionEnd();
	if (endpos > reg2.pos2)
	  reg2.pos2 = endpos;
      }
    
    cluster = toRegionString();
    
  }
  
  void DiscordantCluster::addMateReads(BamReadVector &bav) 
  { 
    
    for (auto& i : bav) {
      std::string sr;
      if (qnames.count(i.Qname())) 
	{
	  std::string tmp = i.GetZTag("SR");
	  i.AddZTag("DC", id);
	  //r_get_Z_tag(i, "SR", tmp);
	  //r_add_Z_tag(i, "DC", id);
	  if (reads.count(tmp) == 0) // only add if this is a mate read
	    mates[tmp] = i;
	}
    }
    
  }
  
  double DiscordantCluster::getMeanMapq(bool mate) const 
  {
    double mean = 0;
    std::vector<int> tmapq;
    if (mate) {
      for (auto& i : mates)
	tmapq.push_back(i.second.MapQuality());
    } else {
      for (auto& i : reads)
	tmapq.push_back(i.second.MapQuality());
    }
    
    if (tmapq.size() > 0)
      mean = accumulate(tmapq.begin(), tmapq.end(), 0.0) / tmapq.size();
    return mean;
  }
  
  double DiscordantCluster::getMeanMapq() const 
  {
    double mean = 0;
    std::vector<int> tmapq;
    for (auto& i : mates)
      tmapq.push_back(i.second.MapQuality());
    for (auto& i : reads)
      tmapq.push_back(i.second.MapQuality());
    
    if (tmapq.size() > 0)
      mean = accumulate(tmapq.begin(), tmapq.end(), 0.0) / tmapq.size();
    return mean;
  }
  
  std::string DiscordantCluster::toRegionString() const 
  {
    //int pos1 = (reg1.strand == '+') ? reg1.pos2 : reg1.pos1;
    //int pos2 = (reg2.strand == '+') ? reg2.pos2 : reg2.pos1;
    int pos1 = reg1.strand ? reg1.pos2 : reg1.pos1;
    int pos2 = reg2.strand ? reg2.pos2 : reg2.pos1;
    
    std::stringstream ss;
    ss << reg1.chr+1 << ":" << pos1 << "(" << reg1.strand << ")" << "-" << 
      reg2.chr+1 << ":" << pos2 << "(" << reg2.strand << ")";
    return ss.str();
    
  }
  
  // define how to print this to stdout
  std::ostream& operator<<(std::ostream& out, const DiscordantCluster& dc) 
  {
    
    out << dc.toRegionString() << " Tcount: " << dc.tcount << 
      " Ncount: "  << dc.ncount << " Mean MAPQ: " 
	<< dc.reads_mapq << " Mean Mate MAPQ: " << dc.mates_mapq;
    /*for (auto& i : dc.reads) {
      std::string tmp;
      i.second->GetTag("SR",tmp);
      out << "   " << i.second->RefID << ":" << i.second->Position << (i.second->IsReverseStrand() ? "(-)" : "(+)") << " - " << i.second->MateRefID << ":" << i.second->MatePosition <<  (i.second->IsMateReverseStrand() ? "(-)" : "(+)") << " " << tmp << endl;
      }*/
    return out;
  }
  
  
  // define how to print to file
  std::string DiscordantCluster::toFileString(bool with_read_names /* false */) const 
  { 
    
    std::string sep = "\t";
    
    // add the reads names (currently off)
    std::string reads_string = "x";
    if (with_read_names) 
      {
      for (auto& i : reads) 
	{
	  std::string tmp = i.second.GetZTag("SR");
	  //r_get_Z_tag(i.second, "SR", tmp);
	  
	  //i.second->GetTag("SR", tmp);
	  reads_string += tmp + ",";
	}
      
      //debug
      if (reads_string.length() == 0)
	reads_string = "x";
      else
	reads_string = reads_string.substr(0,reads_string.length() - 1);
      }
    
    std::stringstream out;
    out << reg1.chr+1 << sep << reg1.pos1 << sep << reg1.strand << sep 
	<< reg2.chr+1 << sep << reg2.pos1 << sep << reg2.strand << sep 
	<< tcount << sep << ncount << sep << reads_mapq << sep 
	<< mates_mapq << sep << reads_string;
    
    return (out.str());
    
  }
  
  // define how to sort theses
  bool DiscordantCluster::operator<(const DiscordantCluster &b) const 
  {
    if (reg1.chr < b.reg1.chr)
      return true;
    if (reg1.pos1 < b.reg1.pos1)
      return true;
    return false;
    
  }
  
}
