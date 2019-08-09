create table goa2(v1 TEXT,v2 TEXT,v3 TEXT,v4 TEXT,v5 TEXT,v6 TEXT,v7 TEXT,v8 TEXT,v9 TEXT);
.import goa_1000.gaf goa;
drop table goa;

create table gene_swiss(v1 TEXT,v2 TEXT);
.import genomeNum_swissprot.id gene_swiss
INSERT INTO gene_swiss (v1, v2) VALUES ("A0A000","icey");


SELECT
  *
FROM
  gene_swiss
LEFT JOIN goa2 ON
  gene_swiss.v2 = goa2.v2;


.separator "\t"
create table gene_swiss_all(v1 TEXT,v2 TEXT,v3 TEXT,v4 TEXT,v5 TEXT,v6 TEXT,v7 TEXT,v8 TEXT,v9 TEXT);

.import goa_uniprot_all.gaf gene_swiss_all

SELECT
  *
FROM
  gene_swiss
LEFT JOIN gene_swiss_all ON
  gene_swiss.v2 = gene_swiss_all.v2;