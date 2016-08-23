#ifndef SEQLIB_CONTIG_PLOT_H__
#define SEQLIB_CONTIG_PLOT_H__

#include "SeqLib/BamRecord.h"

namespace SeqLib {
  
  /** Object for creating ASCII alignment plots
   */
  class SeqPlot {

  public:
    
    /** Create an empty plot */
    SeqPlot() {}

    /** Plot read by stacking them in an IGV-like view */
    std::string PlotAlignmentRecords(const BamRecordVector& brv) const;

    /** Set the view */
    void SetView(const GenomicRegion& g) { m_view = g; }

    /** Set the padding (default is 5) */
    void SetPadding(int p) { m_pad = p; }

  private: 

    // reads that align to the contig
    BamRecordVector m_reads;

    // view window
    GenomicRegion m_view;

    // padding when placing reads
    int m_pad = 5;

  };

  /** A single plotted read */
struct PlottedRead {

  int pos;
  std::string seq;
  std::string info;
  
  bool operator<(const PlottedRead& pr) const {
    return (pos < pr.pos);
  }

};

typedef std::vector<PlottedRead> PlottedReadVector;

/** A plotted line */
struct PlottedReadLine {

  std::vector<PlottedRead*> read_vec;
  int available = 0;
  int contig_len = 0;
  
  int pad = 5;

  void addRead(PlottedRead *r) {
    read_vec.push_back(r);
    available = r->pos + r->seq.length() + pad;
  }

  bool readFits(PlottedRead &r) {
    return (r.pos >= available);
  }

  friend std::ostream& operator<<(std::ostream& out, const PlottedReadLine &r) {
    int last_loc = 0;
    for (auto& i : r.read_vec) {
      assert(i->pos - last_loc >= 0);
      out << std::string(i->pos - last_loc, ' ') << i->seq;
      last_loc = i->pos + i->seq.length();
    }
    int name_buff = r.contig_len - last_loc;
    assert(name_buff < 1e6);
    out << std::string(std::max(name_buff, 5), ' ');
    for (auto& i : r.read_vec) { // add the data
      out << i->info << ",";
    }
    return out;
  }

};

typedef std::vector<PlottedReadLine> PlottedReadLineVector;


}



#endif
