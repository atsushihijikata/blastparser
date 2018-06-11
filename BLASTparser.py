'''
BLASTparser
author: atsushi@nagahama-i-bio.ac.jp
date: Dec. 15, 2017
update: Jun. 11, 2018
'''

class Subject(object):
    '''
    Attributes for BLAST output with the tab-separated format
    '''
    subject_attrs = ('query', 'subject', 'identity', 'match', 'gap', 'mismatch',
                     'qst', 'qen', 'tst', 'ten', 'evalue', 'score')
    subject_types = ('str', 'str', 'float', 'int', 'int', 'int', 'int', 'int',
                     'int', 'int', 'float', 'float')

    def __init__(self, hit):
        for attr, atype, value in zip(self.subject_attrs, self.subject_types, hit):
            if atype == "int":
                value = int(value)
            elif atype == "float":
                value = float(value)
            setattr(self, attr, value)
        self.is_nr = True
        if self.query.startswith('ref|'):
            self.query = self.query.split('|')[1]
        if self.subject.startswith('ref|'):
            self.subject = self.subject.split('|')[1]

IDEN = 30.0
class BLASTparse(object):
    def __init__(self, blastout, nr=False, ordered=True, identity=IDEN):
        self.blastout = blastout
        self.hit_list = []
        self.nr = nr
        self.is_ordered = ordered
        self.threshold = identity

    def parse(self):
        for bo in self.blastout:
            if float(bo[2]) < self.threshold:
                continue
            self.hit_list.append(Subject(bo))
        self._find_split()
        self._find_overlap()
        if self.nr:
            self._find_representative()
        if self.is_ordered:
            self.hit_list = sorted(self.hit_list, key=lambda x:x.qst)

    def _is_split_subject(self, hit1, hit2):
        '''
        query: =========================
        hit1:  ============
        hit2:              ============.
        '''
        #if (hit1.tst < hit2.tst and hit1.ten < hit2.tst) or\
           #(hit1.tst > hit2.tst and hit1.tst > hit2.ten):
        if hit1.ten-5 <= hit2.tst or hit2.ten-5 <= hit1.tst:
            return True
        return False

    def _is_overlap(self, hit1, hit2):
        if hit1.qen >= hit2.qst and hit1.qst <= hit2.qen: # overlap
            if hit2.qen-hit2.qst < 0:
                return False
            if (hit1.qen-hit1.qst+1)/(hit2.qen-hit2.qst+1) <= 1.1 or\
               (hit1.qen-hit1.qst+1)/(hit2.qen-hit2.qst+1) >= 0.9:
                   return True
        else:
            return False

    def _find_overlap(self):
        for i, hit1 in enumerate(self.hit_list):
            for j, hit2 in enumerate(self.hit_list):
                if i == j or hit1.subject != hit2.subject:
                    continue
                if self._is_overlap(hit1, hit2):
                    if hit1.identity >= hit2.identity:
                        hit2.is_nr = False
                    else:
                        hit1.is_nr = False

    def _merge_hits(self, hit1, hit2):
        if int(hit1.ten) < int(hit2.ten):
            hit1.ten = hit2.ten
            hit1.qen = hit2.qen
        else:
            hit1.tst = hit2.tst
            hit1.qst = hit2.qst
        hit2.is_nr = False

    def _find_split(self):
        '''
        The query hits multiple times to a single subject, it means two possibilities:
        1. The query sequence contains multi copies of the subject,
        2. The subject has a long insertion (ex. in case of a fusion protein)
        In the second case, the alignment of query and subject should be done globally.
        '''
        for i, hit1 in enumerate(self.hit_list):
            for j, hit2 in enumerate(self.hit_list):
                if i >= j or hit1.subject != hit2.subject:
                    continue
                if self._is_split_subject(hit1, hit2):
                    self._merge_hits(hit1, hit2)
                #print(hit1.subject, hit2.subject)

    def _find_representative(self):
        for i, hit1 in enumerate(self.hit_list):
            if not hit1.is_nr:
                continue
            for j, hit2 in enumerate(self.hit_list):
                if i >= j or not hit2.is_nr:
                    continue
                if hit1.qst > hit2.qen or hit1.qen < hit2.qst: # if True, no-overlap
                    continue
                cs, ce = 0, 0
                if hit1.qst >= hit2.qst:
                    cs = hit2.qst
                else:
                    cs = hit1.qst
                if hit1.qen >= hit2.qen:
                    ce = hit2.qen
                else:
                    ce = hit1.qen
                if hit1.qst < hit2.qst and hit1.qen > hit2.qen:
                    '''
                    hit1: =============
                    hit2:   =======
                    '''
                    hit2.is_nr = False

                if (hit2.qen-hit2.qst+1) == 0:
                    continue

                if float(abs(ce-cs)+1)/(hit2.qen-hit2.qst+1) > 0.8:
                    '''
                    hit1: =======......
                    hit2: ....###======
                      ###/(###=====) > 0.8
                    '''
                    hit2.is_nr = False

if __name__ == '__main__':
    import sys
    blast_file = sys.argv[1]
    blastout = [line.rstrip().split('\t') for line in open(blast_file)]

    obj = BLASTparse(blastout, nr=True)

    obj.parse()

    for d in obj.hit_list:
        if d.is_nr:
            print("\t".join(map(str, [d.query, d.subject, d.is_nr, d.qst, d.qen, d.tst, d.ten, d.identity])))
