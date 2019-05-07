CREATE TABLE hmmer_results (
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
    descr text
);

ALTER TABLE asmash.hmmer_results OWNER TO biosql;

ALTER TABLE ONLY asmash.hmmer_results
    ADD CONSTRAINT hmmer_results_pkey PRIMARY KEY (target_name, query_name, hmm_from, hmm_to, env_from, env_to, domain_evalue, domain_score, domain_bias);

CREATE INDEX hmmer_results_index ON asmash.hmmer_results USING btree (query_name, acc, fullevalue);

CREATE INDEX hmmer_results_protid ON asmash.hmmer_results USING btree (acc);

