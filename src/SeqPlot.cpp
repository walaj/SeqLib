#include "SeqLib/SeqPlot.h"

namespace SeqLib {

std::string SeqPlot::PlotAlignmentRecords(const BamRecordVector& brv) const {

  PlottedReadVector plot_vec;

  for (auto& i : brv) {
    
    // get the position in the view window
    if (i.ChrID() != m_view.chr)
      continue;

    int pos = i.Position() - m_view.pos1;
    if (pos < 0)
      continue;

    if (i.PositionEnd() > m_view.pos2)
      continue;

    // plot with gaps
    std::string tseq = i.Sequence();
    std::string gapped_seq;

    size_t p = i.AlignmentPosition(); // move along on sequence, starting at first non-clipped base
    for (auto& c : i.GetCigar()) {
      if (c.Type() == 'M') { // 
	assert(p + c.Length() <= tseq.length());
	gapped_seq += tseq.substr(p, c.Length());
      } else if (c.Type() == 'D') {
	gapped_seq += std::string(c.Length(), '-');
      }

      if (c.Type() == 'I' || c.Type() == 'M')
	p += c.Length();
    }

    std::stringstream msg;
    msg << i.Qname() << ">>>" << (i.ChrID() + 1) << ":" << i.Position();
      
    // add to the read plot
    plot_vec.push_back({pos, gapped_seq, msg.str()});
    
  }

  // sort them
  std::sort(plot_vec.begin(), plot_vec.end());

  // make a list of lines
  PlottedReadLineVector line_vec;

  // plot the reads from the ReadPlot vector
  for (auto& i : plot_vec) {
    bool found = false;
    for (auto& j : line_vec) {
      if (j.readFits(i)) { // it fits here
	j.addRead(&i);
	found = true;
	break;
      }
    }
    if (!found) { // didn't fit anywhere, so make a new line
      PlottedReadLine prl;
      prl.pad = m_pad;
      prl.contig_len = m_view.width(); //ac.getSequence().length();
      prl.addRead(&i);
      line_vec.push_back(prl);
    }
  }

  std::stringstream ss;
  
  // plot the lines. Add contig identifier to each
  for (auto& i : line_vec) 
    ss << i << std::endl;

  return ss.str();

  
}


}
