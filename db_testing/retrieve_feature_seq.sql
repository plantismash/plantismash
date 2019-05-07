SELECT t1.seqfeature_id,t1.bioentry_id,t2.start_pos, t2.end_pos, t2.strand, t4.value note, substring(t6.seq, t2.start_pos,t2.end_pos-t2.start_pos+1) seq
FROM seqfeature t1 inner join location t2
on t1.seqfeature_id=t2.seqfeature_id 
inner join term t3 on t1.type_term_id=t3.term_id
inner join seqfeature_qualifier_value t4
on t1.seqfeature_id=t4.seqfeature_id
inner join term t5 on t4.term_id=t5.term_id
inner join biosequence t6 on t1.bioentry_id=t6.bioentry_id
where t3.name='cluster' and t5.name='note' 
limit 10