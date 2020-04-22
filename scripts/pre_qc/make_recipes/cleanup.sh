set -euxo pipefail

		
trash=$1
mkdir -p ${trash}
mvDirToTrash() {
	local dir=$1
	local dest=$2
	if [ -d ${dir} ]; then
		echo -e "moving ${dir} to ${dest}"
		mv ${dir} ${dest}
	fi
}
mvDirToTrash ${data_derived} ${trash}
mvDirToTrash ${data_raw} ${trash}
mvDirToTrash ${clean} ${trash}
mvDirToTrash ${prelim} ${trash}
mvDirToTrash ${logs} ${trash}
mvDirToTrash ${scripts} ${trash}
mvRawDataBack() {
	local file=$1
	local dest=$2
	if [ -f "${file}.bed" ]; then
		echo "moving ${file} to ${dest}"
		mkdir -p ${dest}
		mv ${file}.bed ${dest}
		mv ${file}.bim ${dest}
		mv ${file}.fam ${dest}
		mv "${file}.log" ${dest}
	fi
}
mvRawDataBack "${trash}/data_derived/mega_main_raw/mega_main0" ${data_derived}/mega_main_raw/
mvRawDataBack "${trash}/data_derived/mega_jhu_raw/mega_jhu0" ${data_derived}/mega_jhu_raw/
mvRawDataBack "${trash}/data_derived/mega_cbr_raw/mega_cbr0" ${data_derived}/mega_cbr_raw/
mvRawDataBack "${trash}/data_derived/mega_additions_raw/mega_additions0" ${data_derived}/mega_additions_raw/ 
