--
-- PostgreSQL database dump
--

-- Dumped from database version 9.3.3
-- Dumped by pg_dump version 9.3.3
-- Started on 2014-05-25 18:28:07 CEST

SET statement_timeout = 0;
SET lock_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SET check_function_bodies = false;
SET client_min_messages = warning;

--
-- TOC entry 7 (class 2615 OID 41531)
-- Name: asmash; Type: SCHEMA; Schema: -; Owner: biosql
--

CREATE SCHEMA asmash;


ALTER SCHEMA asmash OWNER TO biosql;

--
-- TOC entry 2442 (class 0 OID 0)
-- Dependencies: 7
-- Name: SCHEMA asmash; Type: COMMENT; Schema: -; Owner: biosql
--

COMMENT ON SCHEMA asmash IS 'Storage of HMMer against nr';


SET search_path = asmash, pg_catalog;

SET default_tablespace = '';

SET default_with_oids = false;

--
-- TOC entry 231 (class 1259 OID 49411)
-- Name: gbacc_gi_mapping; Type: TABLE; Schema: asmash; Owner: biosql; Tablespace: 
--

CREATE TABLE gbacc_gi_mapping (
    genbank_id text,
    gi text NOT NULL
);


ALTER TABLE asmash.gbacc_gi_mapping OWNER TO biosql;

--
-- TOC entry 230 (class 1259 OID 49403)
-- Name: gbacc_protid_mapping; Type: TABLE; Schema: asmash; Owner: biosql; Tablespace: 
--

CREATE TABLE gbacc_protid_mapping (
    genbank_id text,
    protein_id text NOT NULL
);


ALTER TABLE asmash.gbacc_protid_mapping OWNER TO biosql;

--
-- TOC entry 229 (class 1259 OID 41532)
-- Name: hmmer_results; Type: TABLE; Schema: asmash; Owner: biosql; Tablespace: 
--

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

--
-- TOC entry 2280 (class 2606 OID 49418)
-- Name: gi; Type: CONSTRAINT; Schema: asmash; Owner: biosql; Tablespace: 
--

ALTER TABLE ONLY gbacc_gi_mapping
    ADD CONSTRAINT gi PRIMARY KEY (gi);


--
-- TOC entry 2282 (class 2606 OID 49422)
-- Name: gi_unique; Type: CONSTRAINT; Schema: asmash; Owner: biosql; Tablespace: 
--

ALTER TABLE ONLY gbacc_gi_mapping
    ADD CONSTRAINT gi_unique UNIQUE (gi);


--
-- TOC entry 2276 (class 2606 OID 41539)
-- Name: hmmer_results_pkey; Type: CONSTRAINT; Schema: asmash; Owner: biosql; Tablespace: 
--

ALTER TABLE ONLY hmmer_results
    ADD CONSTRAINT hmmer_results_pkey PRIMARY KEY (target_name, query_name, hmm_from, hmm_to, env_from, env_to, domain_evalue, domain_score, domain_bias);


--
-- TOC entry 2278 (class 2606 OID 49410)
-- Name: protein_id; Type: CONSTRAINT; Schema: asmash; Owner: biosql; Tablespace: 
--

ALTER TABLE ONLY gbacc_protid_mapping
    ADD CONSTRAINT protein_id PRIMARY KEY (protein_id);


--
-- TOC entry 2274 (class 1259 OID 41540)
-- Name: hmmer_results_index; Type: INDEX; Schema: asmash; Owner: biosql; Tablespace: 
--

CREATE INDEX hmmer_results_index ON hmmer_results USING btree (query_name, acc, fullevalue);


--
-- TOC entry 2443 (class 0 OID 0)
-- Dependencies: 229
-- Name: hmmer_results; Type: ACL; Schema: asmash; Owner: biosql
--

REVOKE ALL ON TABLE hmmer_results FROM PUBLIC;
REVOKE ALL ON TABLE hmmer_results FROM biosql;
GRANT ALL ON TABLE hmmer_results TO biosql;
GRANT ALL ON TABLE hmmer_results TO PUBLIC;


-- Completed on 2014-05-25 18:28:07 CEST

--
-- PostgreSQL database dump complete
--

