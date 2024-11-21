PlantiSMASH - Plant Secondary Metabolites Analysis Shell
===================================================================

# 1. Installing PlantiSMASH 2.0

## 1.1. Create a Linux Environment in Windows
Skip this step if using a Linux server. Refer to the [WSL setup guide](https://medium.com/@larysha.rothmann/getting-started-with-wsl-for-bioinformatics-setting-up-a-fully-functional-and-pretty-53b9c79b5380).

1. Open the Windows command prompt: Press `Win Key`, search for ‘command prompt’, and select ‘Run as administrator’. Click ‘Yes’ when prompted by Windows to allow the app to make changes.
2. Type in the command prompt: `wsl.exe — install`.
3. Wait for the installation to finish, and type ‘Yes’ if prompted.
4. Restart your computer.
5. Click on Ubuntu.exe to register Ubuntu.
6. To have root as the administrator every time you start, use: `sudo -s`.
7. If any issues arise, unregister Ubuntu in PowerShell with: `wsl --unregister Ubuntu`, then click Ubuntu.exe to register Ubuntu again.

## 1.2. Install Conda
Click Ubuntu.exe or visit your Linux server terminal. Refer to the [Conda manual](https://git.wur.nl/bioinformatics/working-at-bif/-/blob/main/04-conda-manual.md?ref_type=heads).

```bash
cd ~            # or your working directory
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -u -p ./miniconda3
./miniconda3/bin/conda init bash
exit
```
After exiting, reopen with Ubuntu.exe and you should see (base).
```bash
source ~/.bashrc
conda config --set auto_activate_base false
source ~/.bashrc
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
# To check
nano ~/.condarc 
which conda 
```

## 1.3. Install PlantiSMASH
```bash
cd /your/work/folder/path    # Replace with your actual path
mamba create -n plantismash python=2.7.15
conda activate plantismash
git clone https://github.com/plantismash/plantismash.git
cd plantismash
python install_dependencies.py   # Takes about 10 minutes
python download_databases.py   # Only need PFAM for fullhmmer
# For usage instructions of plantismash, use the command:
python run_antismash.py -h
```
## 1.4.1. Run PlantiSMASH on a NCBI referent genome (genebank format) to test
```bash
mkdir Arabidopsis_thaliana
cd Arabidopsis_thaliana
datasets download genome accession GCF_000001735.4 --include gbff
unzip ncbi_dataset.zip
python ../run_antismash.py --clusterblast --knownclusterblast --verbose --debug --limit -1 --taxon plants --outputfolder result/ ncbi_dataset/data/GCF_000001735.4/genomic.gbff
# --clusterblast --knownclusterblast are optional
```
## 1.4.2. For using a genome with gff3 + fasta
```bash
python run_antismash.py --verbose --debug --limit -1 --taxon plants --outputfolder result/ --use_phase --gff3 path/to/gff3/file path/to/fasta/file
# Please check the error message. The genome names in the gff3 file may differ from those in the fasta file, causing an error.
```
## Notion
It is recommended to change the output folder or delete it every time you run PlantiSMASH. Because when using the same output folder, files from the previous run may be partially preserved (only rewriting files with the same names).

# 2. Making your rules
This module is in `antismash\generic_modules\hmm_detection`

## 2.1. How rules find BGCs?
### 2.1.1 Identify domains of proteins
HMMs are used to do that by running hmmerscan. The HMMs files, `cluster_rules.txt`, `hmmdetails.txt`, and `filterhmmdetails.txt`. for plant are in `plant`. The `hmmdetails.txt` controls which HMMs will be used (4th column) and the bitscore cutoff (3rd column) to filter hmmerscan results. Usually the bitscore cutoff is `-1` equal to no filter. 

The result recorded in the output .gbk files. To only show the HMMs with highest bitscore from matches on same proteins sequence range, add them in `filterhmmdetails.txt`. For example, the same domain range match HMM UDPGT and UDPGT_2 , but the output will only show the one with highest bitscore.

Note: Another module by using the command `--full-hmmer` will use Pfam-A.hmm to identify any kind of domains. But the results only recorded in the output .gbk files and did not used in other module.

### 2.1.2. Clustering the neighboring genes with identified domains HMMs matches
This involves the `essential` and `Cutoff` of a rule in `cluster_rules.txt`.
A rule is usually formed as follows:
```bash
Name	        Rule	                             Cutoff (in kb)	Extension (in kb)
Product type	minimum(3,[required],[core_list])	5	           1
```

The `core_list` is a list of HMMs names of domains determining product type, such as `Chal _ stic _ synt _ C/Chal _ stic _ synt _ N` for the rule `polyketide`. Here the two HMMs names are joined by `/`, representing the saved cluster contains at least one of them. If joined by `,` , the cluster containing both domains will be saved.

In `required`, each HMMs name is joined by `,` , representing the saved cluster contains at least one of them. The `required` almost is the list of HMMs record in `hmmdetails.txt` because the clustering starts from the gene with a HMM (in the list) match. Then check the HMM matches of neighboring genes whether in the `required`. If a neighboring gene with the match belongs to the `required`, then add this gene and check the neighbors of it and so on until no new gene adding in the cluster. 

The span of left and right “neighborhood” of a gene on the conift or on the chromosome is calculated dynamically by function `get_dynamic_cutoff_multiplier`.

Simply, `left span = right span = Cutoff*(the span of the nearest ten genes)/10`

In idea situation and `Cutoff = 5`, the neighboring genes are the nearest ten genes. But if one of the nearest is way far away then it will be in the ‘neighborhood’. You can increase the cutoff to consider more neighbor genes. When cutoff = 10, then the neighboring genes must include the nearest 10 genes.

### 2.1.3. Filter the clusters
The domains composition of the cluster contains at least two different domains (`--min-domain-number  2`  is default) and meets `core_list`.

The cluster contains at least three genes with `required` HMMs matches (set by the number behind `minimum(`) and they sharing the similarity below 50% (`--cdh-cutoff 0.5` is default). The similarity is calculated by CD-HIT. Maybe meet the error of memory, the default is `--cdh-memory  2000`.

### 2.1.4. output
The cluster also contain genes in the `Extension` of the both ends biosynthetic genes. The cluster type is same to the name of rule found and saved it. 

## 2.2. How to implant a rule to find a cluster with a certain domain composition
### 2.2.1. Check the `hmmdetails.txt` whether containing the HMMs involved in

If not, add HMMs files in `plant` and update `hmmdetails.txt`
### 2.2.2. Make the rule

`required` can only include few HMMs names (even just `core_list`) if want no others between genes with `required` matches. 

If rules respond to same cluster type, can use `or` to link. Check rule `phenolamide` in `cluster_rules.txt` as example.

In `cluster_rules.txt`, the `cyclopeptide` rule will save every cluster contains gene with `BURP` HMMs match. Can follow this format to make rule.

There are other ways to form rules, can check the `__init__.py`

### 2.2.3. Add the rule in `cluster_rules.txt`
Keep the first line; can delete other rules; do not leave blank lines

### 2.2.4. Run plantiSMASH with some commands to change the default (if need)
Such as `--cdh-cutoff`, `--min-domain-number`, `--cdh-memory`

### 2.2.5. Beautify (add color) the output on the web page (optional)

Give violet background to `alkaloid` clusters on overview page
`antismash/output_modules/html/css/style.css`
```bash
.alkaloid {
  background-color: violet;
}
.alkaloid a {
  color: black;
}
````
Set cluster legend of genes with pyridoxal synthase domains HMMs matches.
`antismash/output_modules/html/js/gene_colors.js`
```bash
{ label: "pyridoxal synthase", color: "#4CAF50", members : ["plants/YjeF_N", "plants/Pyridox_oxidase", "plants/PNPOx_C"] },
````

# 3. the Subgroup Module

## 3.1. What is the Subgroup Module?

The subgroup module is located in `antismash/generic_modules/subgroup`. A subgroup is defined as a clade of a family phylogenetic tree, which contains members sharing similar functions. This module enables the subgroup of protein sequences from various families based on the results of `hmmer` scans using domain pHMMs. The families it covers include the Cellulose synthase (CSLs) family (Chung et al., 2020; Jozwiak et al., 2020), UDP-glucuronosyltransferase (UGTs) family (Louveau & Osbourn, 2019), short-chain dehydrogenases/reductases (SDRs) family (Moummou et al., 2012), and oxidosqualene cyclase (OSCs) family (unpublished work).

## 3.2. Configuring the Subgroup Module

You can specify which families to subgroup by altering the ‘Enable’ column from ‘Y’ to ‘N’ in `Subgroup_Model.txt`. Each family folder (`Cellulose_synt`, `OSC`, `SDR`, `UDPGT`) contains HMM profiles of subgroups. Note that for SDR, the HMMs are not based on a phylogenetic tree.

The `Family model name(s)` column in the configuration file designates which domain pHMMs can be associated with a particular family. The `Subgroup Model` and `Match type` columns determine what is displayed on the web pages for each cluster, in the overview, and within .gbk files.

Similar ideas in `Subgroup_Tree.txt`, the reference package (in family folders) containing the referent tree is for the pplacer to place the sequence on the tree.

## Using the Subgroup Module

Commands for the module include:

- `--disable_subgroup` to disable identifying the subgroup.
- `--disable_treesvg` to disable the creation of SVG pictures of subgrouping trees to save time.
- `--subgroup_inputpath SUBGROUP_INPUTPATH` to specify the path to a folder with the same structure as the subgroup folder in `antismash/generic_modules`.

## 3.3. Customizing the Subgroup Module

Take the Cellulose synthase family as an example. If a subgroup of the Cellulose synthase on the family tree has been well studied, you can collect those sequences, divide them into subgroups, and create HMMs using `hmmbuild`. Ensure to add the corresponding information in `Subgroup_Model.txt`.

## Notion
The .hmm file name and the NAME in the .hmm file content should match.

The reference package, which contains the reference tree for `pplacer` to place the sequence on, is created using `taxtastic` with the following components:

- Alignment of reference sequences.    `CSLs_all_to_check_with_add.afa`
- The reference tree based on the alignment. `RAxML_bipartitionsBranchLabels.CSLs_all_to_check_with_add.newick`
- The log file from the tool (iqtree, raxmlng, raxml, fasttree,phyml) used to create the tree. `RAxML_info.CSLs_all_to_check_with_add`
- A table containing node names and their corresponding subgroup names. `node-subgroup.txt`
```bash
cd /your/work/folder/path    # Replace with your actual path
git clone https://github.com/fhcrc/taxtastic.git
cd taxtastic
mamba  create -n plantismash  python=3
conda activate python3
source taxtastic-env/bin/activate
pip install .
cd ../plantismash/antismash/generic_modules/subgroup/Cellulose_synt/Cellulose_synt.refpkg
taxit create -P  Cellulose_synt.refpkg  -l plant_Cellulose_synt  --aln-fasta  CSLs_all_to_check_with_add.afa  --tree-stats  RAxML_info.CSLs_all_to_check_with_add    --tree-file   RAxML_bipartitions.CSLs_all_to_check_with_add  --seq-info node-subgroup.txt
```
`Taxtastic` generates `CONTENTS.json` recording the files and `phylo_model.json` recording information from the log file. Verify `phylo_model.json` to ensure the parameters match those in the log file.  If it is not, it can be changed directly on it.

For example, change `"subs_model": "AUTO"` to `"subs_model": "LG"`.  The same strategy can be used when using a tree-building tool log that is not accepted by `taxtastic`.

After creating the reference package, add the corresponding information to `Subgroup_Tree.txt`.

## Notion
Module uses ete3 to parse the tree file. But to ete3, 0.003847[100] format is not supported. bootstrap values in newick format should look like: ')100:0.003847'.

## reference for the Subgroup Module
Chung, S. Y., Seki, H., Fujisawa, Y., Shimoda, Y., Hiraga, S., Nomura, Y., Saito, K., Ishimoto, M., & Muranaka, T. (2020). A cellulose synthase-derived enzyme catalyses 3-O-glucuronosylation in saponin biosynthesis. Nature Communications 2020 11:1, 11(1), 1–11. https://doi.org/10.1038/s41467-020-19399-0


Jozwiak, A., Sonawane, P. D., Panda, S., Garagounis, C., Papadopoulou, K. K., Abebie, B., Massalha, H., Almekias-Siegl, E., Scherf, T., & Aharoni, A. (2020). Plant terpenoid metabolism co-opts a component of the cell wall biosynthesis machinery. Nature Chemical Biology 2020 16:7, 16(7), 740–748. https://doi.org/10.1038/s41589-020-0541-x


Louveau, T., & Osbourn, A. (2019). The Sweet Side of Plant-Specialized Metabolism. Cold Spring Harbor Perspectives in Biology, 11(12), a034744. https://doi.org/10.1101/CSHPERSPECT.A034744


Moummou, H., Kallberg, Y., Tonfack, L. B., Persson, B., & van der Rest, B. (2012). The Plant Short-Chain Dehydrogenase (SDR) superfamily: Genome-wide inventory and diversification patterns. BMC Plant Biology, 12(1), 1–17. https://doi.org/10.1186/1471-2229-12-219/FIGURES/7

# 4. Update the ClusterBlast Module databas
## 4.1. What is the ClusterBlast Module?
Activate it with `--clusterblast`

It will find the homologous clusters stored in the database

The scripts are in `antismash/generic_modules/clusterblast` also with database files, including:`plantgeneclusters.txt, plantgeneclusterprots.fasta`

Use `--clusterblastdir` to specify the database directory which you want to use.

## 4.2. How to update the ClusterBlast Module database?
1. change the minimum number of each rule in `plant/cluster_rules.txt` to 2 (so will save clusters only with 2 genes);
2. download the genomes you want to use to make the database, for example NCBI reference genomes of Streptophyta:
```bash
datasets download genome taxon 35493 --annotated --reference --include gbff  --dehydrated  --filename streptophyta_ref_anno.zip
unzip streptophyta_ref_anno.zip  -d ncbi_plant_ref_anno
datasets rehydrate --directory ncbi_plant_ref_anno
# check whether get all genomes mentioned in fetch.txt
find ncbi_plant_ref_anno/ncbi_dataset/data -mindepth 1 -maxdepth 1 -type d ! -exec sh -c 'ls "{}" | grep -q "genomic.gbff"' \; -print
```
3. run plantiSMASH with `--clusterblast` to get the clusters, for example for using those Streptophyta genomes:
```bash 
nohup bash -c 'find ncbi_dataset/data -mindepth 1 -maxdepth 1 -type d -exec bash -c "python2 ../plantismash/run_antismash.py  --cpus 90 --cdh-memory 64000 --cdh-cutoff 0.9  --update_clusterblast --verbose  --disable-svg   --disable-html  --disable-xls --disable_specific_modules  --disable_subgroup  --disable-genbank   --limit -1 --taxon plants  --outputfolder {}/for_blast_db  {}/*.*" \; > clusterblast_database_making.log 2>&1' &
```
This will update `plantgeneclusters.txt, plantgeneclusterprots.fasta` in every run

Can add `--clusterblastdir` to specify the directory which you want to save the database.

4. get the numbers of clusters found in each genome:
```bash
grep -E 'INFO: Numbers of clusters' clusterblast_database_making.log > clusters_numbers.txt
```
