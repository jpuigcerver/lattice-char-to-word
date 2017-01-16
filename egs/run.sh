#!/bin/bash
set -e;
export LC_NUMERIC=C;

# Directory where the prepare.sh script is placed.
SDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)";
[ "$SDIR" != "$PWD" ] &&
echo "Please, run this script from the egs directory!" >&2 &&
exit 1;

function get_lattice_nbest () {
  lattice-to-nbest --n=3 ark:- ark:- 2> /dev/null |
  nbest-to-linear ark:- ark:/dev/null ark,t:- ark,t:- 2> /dev/null |
  awk '{
    if (NR % 2 == 0) {
      COST[$1] = $2;
    } else {
      PATH[$1] = $2;
      for(i=3;i<=NF;++i) {
        PATH[$1] = PATH[$1]" "$i;
      }
    }
  }END{
    for (n in PATH) {
      print COST[n], PATH[n];
    }
  }' | sort -g
}

function nbest_word_to_char () {
  awk -v MF="$1" '
BEGIN{
  while((getline < MF) > 0) { M[$2] = $1; }
}{
  for (i=2; i<=NF; ++i)  {
    if ($i in M) {
      $i = M[$i];
    } else {
      print "Symbol "$i" was not found in the table!" > "/dev/stderr";
      exit(1);
    }
  }
  print;
}'
}

# Determine which lattice-info command should be used
lattice_info_cmd="$(which lattice-info 2> /dev/null)" ||
lattice_info_cmd=lattice-info/lattice-info;

if [ ! -f "$lattice_info_cmd" ]; then
  echo "lattice-info was not found in your system, compiling!" >&2;
  [ -d lattice-info/.git ] ||
  git clone https://github.com/jpuigcerver/lattice-info.git || exit 1;
  ( cd lattice-info && make && cd .. ) || exit 1;
fi;

# Convert character-level to word-level lattice
../lattice-char-to-word \
  --save-isymbols=lattice.word.sym --save-osymbols=/dev/null \
  3 ark:lattice.char.txt ark,t:lattice.word.txt;

# Both lattices should have the same number of paths!
echo "CHECK NUMBER OF PATHS...";
"$lattice_info_cmd" ark:lattice.char.txt | grep "avg. of paths";
"$lattice_info_cmd" ark:lattice.word.txt | grep "avg. of paths";
echo "";

# All paths should have the same cost!
echo "CHECK PATHS AND COSTS...";
echo "Paths from char lattice:";
cat lattice.char.txt | get_lattice_nbest;
echo "Paths from word lattice:";
cat lattice.word.txt | get_lattice_nbest | nbest_word_to_char lattice.word.sym;
