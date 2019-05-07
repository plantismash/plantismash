DROP TABLE asmash.clusterblast_table;
CREATE TABLE asmash.clusterblast_table
(
  gi text NOT NULL,
  prot_id text NOT NULL,
  subject_string text NOT NULL,
  perc_id numeric NOT NULL,
  aln_length integer NOT NULL,
  no_mism smallint NOT NULL,
  no_gaps smallint NOT NULL,
  query_start integer NOT NULL,
  query_end integer NOT NULL,
  subject_start integer NOT NULL,
  subject_end integer NOT NULL,
  evalue numeric NOT NULL,
  bitscore numeric NOT NULL,
  CONSTRAINT clusterblast_table_pkey PRIMARY KEY (gi, prot_id, subject_string, query_start, query_end, subject_start, subject_end)
);

--DROP INDEX clusterblast_table_index;
CREATE INDEX clusterblast_table_index
  ON asmash.clusterblast_table
  USING btree
  (gi, prot_id, evalue, bitscore);
  
--DROP INDEX clusterblast_table_protid
  CREATE INDEX clusterblast_table_protid 
  ON asmash.clusterblast_table 
  USING btree 
  (prot_id);