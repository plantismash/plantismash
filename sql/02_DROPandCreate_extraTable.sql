SET statement_timeout = 0;
SET lock_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SET check_function_bodies = false;
SET client_min_messages = warning;



SET search_path = antismashExtradata, pg_catalog;

SET default_tablespace = '';

SET default_with_oids = false;


CREATE TABLE antismashExtradata.extradata (
    genbank_id text,
    secondary_id text,
    asdata text NOT NULL
);

ALTER TABLE extradata OWNER TO biosql;
ALTER TABLE ONLY extradata  ADD CONSTRAINT genbank_id PRIMARY KEY (genbank_id, secondary_id);
ALTER TABLE ONLY extradata  ADD CONSTRAINT genbank_id_unique UNIQUE (genbank_id, secondary_id);
CREATE INDEX extradataIndex ON extradata USING btree (genbank_id);
CREATE INDEX secondarydataIndex ON extradata USING btree (genbank_id, secondary_id);

REVOKE ALL ON TABLE extradata FROM PUBLIC;
REVOKE ALL ON TABLE extradata FROM biosql;
GRANT ALL ON TABLE extradata TO biosql;
