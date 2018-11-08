esearch -db gene -query "deuteranopia" |
  efetch -format xml |
  xtract -pattern Entrezgene \
    -element Gene-track_geneid Gene-ref_locus \
    -sep "|" -element Gene-ref_syn_E



 esearch -db gene -query "XP_015832938.1" |
  efilter -status alive | efetch -format docsum |
  xtract -pattern DocumentSummary -NAME Name -DESC Description \
    -block GenomicInfoType -element "&NAME" "&DESC" 


cat all.list|cut -d" " -f 1|sort -k1,1 -u > gene.list #118


cat gene.list|xargs -t -I {} esearch -db gene -query "{}" |
  efilter -status alive | efetch -format docsum |
  xtract -pattern DocumentSummary -NAME Name -DESC Description \
    -block GenomicInfoType -element "&NAME" "&DESC" >list


cat gene.list | xargs -t -i sh -c 'esearch -db gene -query "{}" |
  efilter -status alive | efetch -format docsum |
  xtract -pattern DocumentSummary -NAME Name -DESC Description \
    -block GenomicInfoType -element "&NAME" "&DESC"' >list









