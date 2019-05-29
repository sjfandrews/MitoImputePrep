#!/usr/bin/env bash
rm in_bash.scripts
rm in_nonbash.scripts
rm in_no.scripts

find . -name "*.sh" -o -name "*.R" -o -name "*.py" -o -name "*.smk" > scripts.list

grep -v "scripts/BASH" scripts.list > scripts_nobash.list
grep "scripts/BASH" scripts.list > scripts_bash.list

echo EXCLUDING files from the BASH folder:
echo
for i in $(cat scripts.list); do
  sed 's/ /\\\ /g' scripts_nobash.list | \
	  xargs grep -l "$(basename $i)" | \
    grep -v $i > scripts.list.tmp
	if [ -s scripts.list.tmp ]; then
    echo $i >> in_nonbash.scripts
	  echo $i is in the following scripts:
	  cat scripts.list.tmp
		rm scripts.list.tmp
    echo
    echo
	fi
done

echo INCLUDING ONLY files from the BASH folder:
echo
for i in $(cat scripts.list); do
  sed 's/ /\\\ /g' scripts_bash.list | \
	  xargs grep -l "$(basename $i)" | \
    grep -v $i > scripts.list.tmp
	if [ -s scripts.list.tmp ]; then
    echo $i >> in_bash.scripts
	  echo $i is in the following scripts:
	  cat scripts.list.tmp
    echo
    echo
	fi
  rm scripts.list.tmp
done

echo ORPHAN SCRIPTS:
echo
for i in $(cat scripts.list); do
  sed 's/ /\\\ /g' scripts.list | \
	  xargs grep -l "$(basename $i)" | \
    grep -v $i > scripts.list.tmp
	if [ ! -s scripts.list.tmp ]; then
    echo $i >> in_no.scripts
	  echo $i
	fi
  rm scripts.list.tmp
done
rm scripts.list
rm scripts_nobash.list
rm scripts_bash.list
