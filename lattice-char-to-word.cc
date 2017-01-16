// MIT License
//
// Copyright (c) 2016 Joan Puigcerver <joapuipe@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include "base/kaldi-common.h"
#include "util/common-utils.h"
#include "fstext/kaldi-fst-io.h"
#include "fstext/fstext-utils.h"
#include "lat/kaldi-lattice.h"
#include "lat/lattice-functions.h"

namespace fst {

template <typename Arc, bool match_input = false>
void ExpandFst(
    const Fst<Arc>& ifst,
    const std::unordered_set<typename Arc::Label>& delim,
    const int32 max_length, MutableFst<Arc>* ofst,
    std::unordered_map<std::vector<typename Arc::Label>,
                       typename Arc::Label>* ilabel_map,
    std::unordered_map<std::vector<typename Arc::Label>,
                       typename Arc::Label>* olabel_map) {
  typedef Fst<Arc> Fst;
  typedef typename Arc::Label Label;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  typedef std::tuple< StateId,
                      StateId,
                      Weight,
                      std::vector<Label>,
                      std::vector<Label> > SState;

  // Reserve symbol for epsilon
  ilabel_map->insert(make_pair(std::vector<Label>{}, 0));
  olabel_map->insert(make_pair(std::vector<Label>{}, 0));

  // Output fst has, at most, as many states as the input fst.
  ofst->DeleteStates();
  for (StateIterator<Fst> siter(ifst); !siter.Done(); siter.Next()) {
    ofst->AddState();
  }

  // Set start state and add it to the stack, since from this state new words
  // may begin.
  std::unordered_set<StateId> add_to_initial_stack;
  ofst->SetStart(ifst.Start());
  add_to_initial_stack.insert(ifst.Start());

  // Traverse all arcs from the fst to detect where new words may start.
  // We will add these states to the stack, and keep the arcs with any of the
  // word delimiters, in order to keep these as words.
  for (StateIterator<Fst> siter(ifst); !siter.Done(); siter.Next()) {
    const StateId s = siter.Value();
    ofst->SetFinal(s, ifst.Final(s));

    for (ArcIterator<Fst> aiter(ifst, s); !aiter.Done(); aiter.Next()) {
      const Arc& arc = aiter.Value();
      const Label& match_label = match_input ? arc.ilabel : arc.olabel;
      // State arc.nextstate is the start of a new word
      if (delim.count(match_label) != 0) {
        // Get ilabel and olabel corresponding to this arc in the output fst
        std::vector<Label> i_l, o_l;
        if (arc.ilabel != 0) i_l.push_back(arc.ilabel);
        if (arc.olabel != 0) o_l.push_back(arc.olabel);
        const Label ilab = ilabel_map->insert(
            make_pair(i_l, ilabel_map->size())).first->second;
        const Label olab = olabel_map->insert(
            make_pair(o_l, olabel_map->size())).first->second;
        // Add arc to the output fst and add it to the stack
        ofst->AddArc(s, Arc(ilab, olab, arc.weight, arc.nextstate));
        add_to_initial_stack.insert(arc.nextstate);
      }
    }
  }

  // Initialize stack
  std::stack<SState> S;
  for (StateId s : add_to_initial_stack) {
    S.push(SState(s, s, Weight::One(),
                  std::vector<Label>{}, std::vector<Label>{}));
  }

  while (!S.empty()) {
    const SState s = S.top(); S.pop();
    const StateId& q0 = std::get<0>(s), q1 = std::get<1>(s);
    const Weight& w = std::get<2>(s);
    const std::vector<Label> ilbls = std::get<3>(s), olbls = std::get<4>(s);

    bool has_delim_arc = false;
    for (ArcIterator<Fst> aiter(ifst, q1); !aiter.Done(); aiter.Next()) {
      const Arc& arc = aiter.Value();
      const Label match_label = match_input ? arc.ilabel : arc.olabel;
      if (delim.count(match_label) == 0) {
        const size_t cur_length = match_input ? ilbls.size() : olbls.size();
        const size_t new_length = cur_length + (match_label == 0 ? 0 : 1);
        if (new_length <= max_length) {
          // Continue expanding from (q0, arc.nextstate)
          // Accumulate weight and input/output labels.
          std::vector<Label> i_l(ilbls), o_l(olbls);
          if (arc.ilabel != 0) i_l.push_back(arc.ilabel);
          if (arc.olabel != 0) o_l.push_back(arc.olabel);
          S.push(SState(q0, arc.nextstate, Times(w, arc.weight), i_l, o_l));
        }
      } else {
        has_delim_arc = true;
      }
    }

    // If q1 is final or has any output arc with a delimiter symbol,
    // add arc from q0 to q1 to the output fst, with the cumulative weight
    // and input/output labels.
    if (q0 != q1 && (has_delim_arc || ifst.Final(q1) != Weight::Zero())) {
      const Label ilab = ilabel_map->insert(
          make_pair(ilbls, ilabel_map->size())).first->second;
      const Label olab = olabel_map->insert(
          make_pair(olbls, olabel_map->size())).first->second;
      ofst->AddArc(q0, Arc(ilab, olab, w, q1));
    }
  }

  // Trim output fst
  Connect(ofst);
}

template <typename M>
void MapToSymbolsTable(
    const M& smap, SymbolTable* stable) {
  typedef typename M::key_type key_type;        // i.e. vector<int>
  typedef typename M::mapped_type mapped_type;  // i.e. int
  for (const auto& k_v : smap) {
    const key_type& k = k_v.first;
    const mapped_type& v = k_v.second;
    std::ostringstream ss;
    for (const mapped_type& c : k) {
      ss << c << "_";
    }
    std::string str = ss.str();
    if (str.empty()) { str = "0"; }  // epsilon
    else { str.pop_back(); }         // pop last _
    KALDI_ASSERT(stable->AddSymbol(str, v) == v);
  }
}

}  // namespace fst

namespace std {

template <>
struct hash<std::vector<kaldi::CompactLatticeArc::Label>> {
  typedef std::vector<kaldi::CompactLatticeArc::Label> argument_type;
  size_t operator()(const argument_type& v) const {
    hash<kaldi::CompactLatticeArc::Label> hashe;
    size_t h = 0;
    for (const auto& s : v) {
      h ^= hashe(s) + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
  }
};

}  // namespace std

int main(int argc, char** argv) {
  try {
    using namespace kaldi;

    const char* usage =
        "Convert character-level lattices into word-level lattices by "
        "expanding the subpaths in between any of two separator symbols.\n"
        "\n"
        "Keep in mind that this lattice expansion has an exponential cost. "
        "For instance, if the set of separator symbols was empty, all paths "
        "from the input lattice would be expanded, so that each arc in the "
        "output lattice would be a full path from the input lattice.\n"
        "\n"
        "Usage: lattice-char-to-word [options] separator-symbols "
        "lat-rspecifier lat-wspecifier\n"
        " e.g.: lattice-char-to-word \"3 4\" ark:1.lat ark:1-words.lat\n";

    ParseOptions po(usage);
    BaseFloat acoustic_scale = 1.0;
    BaseFloat graph_scale = 1.0;
    BaseFloat beam = std::numeric_limits<BaseFloat>::infinity();
    int32 max_length = std::numeric_limits<int32>::max();
    std::string save_isymbols = "";
    std::string save_osymbols = "";

    po.Register("acoustic-scale", &acoustic_scale,
                "Scaling factor for acoustic likelihoods in the lattices.");
    po.Register("graph-scale", &graph_scale,
                "Scaling factor for graph probabilities in the lattices.");
    po.Register("beam", &beam, "Pruning beam (applied after acoustic scaling "
                "and adding the insertion penalty).");
    po.Register("save-isymbols", &save_isymbols, "");
    po.Register("save-osymbols", &save_osymbols, "");
    po.Read(argc, argv);

    if (po.NumArgs() != 3) {
      po.PrintUsage();
      exit(1);
    }

    // Parse delimiter symbols from arguments
    typedef CompactLatticeArc::Label Label;
    std::unordered_set<Label> delimiter_symbols;
    {
      std::istringstream delimiter_symbols_iss(po.GetArg(1));
      Label tmp;
      while (delimiter_symbols_iss >> tmp) {
        if (tmp == 0) {
          KALDI_ERR << "Epsilon (0) cannot be a delimiter symbol!";
        }
        delimiter_symbols.insert(tmp);
      }
    }

    if (graph_scale <= 0.0 || acoustic_scale <= 0.0) {
      KALDI_ERR << "--acoustic-scale and --graph-scale must be strictly greater than 0.0!";
    }
    // Scaling scores
    std::vector<std::vector<double> > scale(2, std::vector<double>{0.0, 0.0});
    scale[0][0] = graph_scale;
    scale[1][1] = acoustic_scale;
    // Inv. scaling scores, used to rescale the output lattices before writing
    std::vector<std::vector<double> > iscale(2, std::vector<double>{0.0, 0.0});
    scale[0][0] = 1.0 / graph_scale;
    scale[1][1] = 1.0 / acoustic_scale;

    const std::string lattice_in_str = po.GetArg(2);
    const std::string lattice_out_str = po.GetArg(3);

    std::unordered_map<std::vector<Label>, Label> ilabel_map;
    std::unordered_map<std::vector<Label>, Label> olabel_map;
    SequentialCompactLatticeReader lattice_reader(lattice_in_str);
    CompactLatticeWriter lattice_writer(lattice_out_str);

    for (; !lattice_reader.Done(); lattice_reader.Next()) {
      const std::string lattice_key = lattice_reader.Key();

      CompactLattice lat = lattice_reader.Value();
      lattice_reader.FreeCurrent();

      // Acoustic scale
      if (acoustic_scale != 1.0 || graph_scale != 1.0)
        fst::ScaleLattice(scale, &lat);

      // Prune lattice
      if (beam != std::numeric_limits<BaseFloat>::infinity())
        PruneLattice(beam, &lat);

      // Put weights into the original scale
      if (acoustic_scale != 1.0 || graph_scale != 1.0)
        fst::ScaleLattice(iscale, &lat);

      // Expand fst to obtain the word-level arcs
      CompactLattice olat;
      if (save_isymbols == "") ilabel_map.clear();
      if (save_osymbols == "") olabel_map.clear();
      fst::ExpandFst(lat, delimiter_symbols, max_length, &olat,
                     &ilabel_map, &olabel_map);

      // Write symbol tables in the output lattice, if no filenames where given
      if (save_isymbols == "") {
        fst::SymbolTable stable;
        MapToSymbolsTable(ilabel_map, &stable);
        olat.SetInputSymbols(&stable);
      }
      if (save_osymbols == "") {
        fst::SymbolTable stable;
        MapToSymbolsTable(olabel_map, &stable);
        olat.SetOutputSymbols(&stable);
      }

      // Write output lattice
      lattice_writer.Write(lattice_key, olat);
    }

    // Save input and output symbol tables into separate files, if these
    // filenames where given
    if (save_isymbols != "") {
      fst::SymbolTable stable;
      MapToSymbolsTable(ilabel_map, &stable);
      KALDI_ASSERT(stable.WriteText(save_isymbols));
    }
    if (save_osymbols != "") {
      fst::SymbolTable stable;
      MapToSymbolsTable(olabel_map, &stable);
      KALDI_ASSERT(stable.WriteText(save_osymbols));
    }

    return 0;
  } catch (const std::exception& e) {
    std::cerr << e.what();
    return 1;
  }
}
