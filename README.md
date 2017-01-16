# lattice-char-to-word

This tool convert character-level lattices into word-level lattices by
expanding the subpaths in between any of two separator symbols.

Keep in mind that this lattice expansion has an exponential cost. For instance,
if the set of separator symbols was empty, all paths from the input lattice
would be expanded, so that each arc in the output lattice would be a full path
from the input lattice.

However, the exponential growth is constrained by the use of separator symbols,
and make the tool practical in real scenarios.

In addition, there are two prunning mechanisms to prevent the output lattices
from exploding:

1. You can prune the character lattices before expanding with the --beam option.
2. You can set a maximum length for the output words with the --max-length
   option. Any path with a word longer than this number of characters will
   be removed from the output path.

```bash
Usage: lattice-char-to-word [options] separator-symbols lat-rspecifier lat-wspecifier
 e.g.: lattice-char-to-word "3 4" ark:1.lat ark:1-words.lat

Options:
  --acoustic-scale            : Scaling factor for acoustic likelihoods in the lattices. (float, default = 1)
  --beam                      : Pruning beam (applied after acoustic scaling and adding the insertion penalty). (float, default = inf)
  --graph-scale               : Scaling factor for graph probabilities in the lattices. (float, default = 1)
  --max-length                : Max. length (in characters) for a word. (int, default = 2147483647)
  --save-symbols              : If given, all lattices will use the same symbol table which will be written to this destination file. If not provided, each lattice contains its own symbol table. (string, default = "")

```
