#!/bin/bash
set -euxo pipefail

group=$1
password=$2
cd ${group}

unzipResult () {

		local i=$1
		if [ -f chr_${i}.zip ]; then
			echo -e "UNZIP CHR${i}"
			7za e -p"${password}" chr_${i}.zip
			wait
			#rm -f chr_${i}.zip
		fi

}


for chr in {1..22}; do
	unzipResult "${chr}"
done
wait

echo -e "FINISHED UNZIPPING ${group}"
