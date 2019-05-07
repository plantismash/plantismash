DROP TABLE asmash.pfam_table;
CREATE TABLE asmash.pfam_table
(
  target_name text NOT NULL,
  acc text NOT NULL,
  tlen integer,
  query_name text NOT NULL,
  acc2 text,
  qlen smallint,
  fullevalue numeric,
  full_score double precision,
  full_bias double precision,
  domain_num smallint,
  domain_of smallint,
  domain_evalue numeric NOT NULL,
  domain_ievalue numeric,
  domain_score double precision NOT NULL,
  domain_bias double precision NOT NULL,
  hmm_from integer NOT NULL,
  hmm_to integer NOT NULL,
  ali_from integer,
  ali_to integer,
  env_from integer NOT NULL,
  env_to integer NOT NULL,
  acc3 double precision,
  descr text,
  CONSTRAINT pfam_table_pkey PRIMARY KEY (target_name, query_name, hmm_from, hmm_to, env_from, env_to, domain_evalue, domain_score, domain_bias)
);

ALTER TABLE asmash.pfam_table
  OWNER TO biosql;

-- DROP INDEX asmash.pfam_table_index;

CREATE INDEX pfam_table_index
  ON asmash.pfam_table
  USING btree
  (query_name, acc, fullevalue);

CREATE INDEX pfam_table_protid 
  ON asmash.pfam_table 
  USING btree 
  (acc);
