SELECT 
  biosequence.length, 
  location.start_pos, 
  location.end_pos, 
  location.strand,
  substring(biosequence.seq,location.start_pos, location.end_pos-location.start_pos+1) seq,
  reverse(compl(substring(biosequence.seq,location.start_pos, location.end_pos-location.start_pos+1))) revseq,
  bioentry.bioentry_id, 
  bioentry.name, 
  bioentry.accession, 
  seqfeature_qualifier_value.value, 
  term.name
FROM 
  public.bioentry, 
  public.biosequence, 
  public.location, 
  public.seqfeature, 
  public.seqfeature_qualifier_value, 
  public.term
WHERE 
  bioentry.bioentry_id = biosequence.bioentry_id AND
  bioentry.bioentry_id = seqfeature.bioentry_id AND
  seqfeature.seqfeature_id = location.seqfeature_id AND
  seqfeature.seqfeature_id = seqfeature_qualifier_value.seqfeature_id AND
  seqfeature_qualifier_value.term_id = term.term_id AND
  bioentry.name = 'Y16952' AND
  seqfeature_qualifier_value.value LIKE 'AMP-binding';
