#!/bin/bash
set -o nounset
set -o errexit
SRCDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

$SRCDIR/../extract-dihedrals.py $SRCDIR/example-input.gro -ignoreres SOL NA CL
