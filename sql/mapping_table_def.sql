SET statement_timeout = 0;
SET lock_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SET check_function_bodies = false;
SET client_min_messages = warning;

SET search_path = asmash, pg_catalog;

SET default_tablespace = '';

SET default_with_oids = false;

CREATE TABLE gbacc_gi_mapping (
    genbank_id text NOT NULL,
    gi text NOT NULL,
    start bigint NOT NULL,
    "end" bigint NOT NULL,
    strand integer NOT NULL
);


ALTER TABLE asmash.gbacc_gi_mapping OWNER TO biosql;

CREATE TABLE gbacc_protid_mapping (
    genbank_id text NOT NULL,
    protein_id text NOT NULL,
    start bigint NOT NULL,
    "end" bigint NOT NULL,
    strand integer NOT NULL
);

ALTER TABLE asmash.gbacc_protid_mapping OWNER TO biosql;

ALTER TABLE ONLY gbacc_gi_mapping
    ADD CONSTRAINT "gi-mapping" PRIMARY KEY (genbank_id, gi);

ALTER TABLE ONLY gbacc_gi_mapping
    ADD CONSTRAINT "gi-unique" UNIQUE (genbank_id, gi);

ALTER TABLE ONLY gbacc_protid_mapping
    ADD CONSTRAINT "prot_id-mapping" PRIMARY KEY (genbank_id, protein_id);

ALTER TABLE ONLY gbacc_protid_mapping
    ADD CONSTRAINT "prot_id-unique" UNIQUE (genbank_id, protein_id);

CREATE INDEX gbacc_protid_mapping_genbank_id 
	ON asmash.gbacc_protid_mapping 
	USING btree
	(genbank_id);
