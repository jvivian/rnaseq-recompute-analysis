class GTFRecord(object):
    def __init__(self, line):
        """
        Converts a GTF record into an python object

        :param line str|list: GTF record
        """
        if isinstance(line, str):
            line = line.strip().split('\t')
        if len(line) != 9:
            msg = '\n'.join(self.attributes)
            raise ValueError('GTF file format has changed. Must have these attributes: {}'.format(msg))
        for (attr, func), value in zip(self.attributes, line):
            setattr(self, attr, func(value))
        for feature in self.attribute.split(';'):
            try:
                attr, value = feature.split()
                setattr(self, attr, value.replace('"', ''))
            except ValueError:
                pass

    attributes = [('seqname', str), ('source', str), ('feature', str), ('start', int), ('end', int),
                  ('score', str), ('strand', str), ('frame', str), ('attribute', str)]

    def __repr__(self):
        if self.feature == 'gene':
            return 'GTF({}, gene:{}, start:{}, length:{})'.format(self.seqname, self.gene_name,
                                                                  self.start, self.end-self.start+1)
        elif self.feature == 'start_codon':
            return 'GTF({}, start codon:{}, start:{}, length:{})'.format(self.seqname,
                                                                         self.gene_name,
                                                                         self.start,
                                                                         self.end - self.start + 1)
        elif self.feature == 'transcript':
            return 'GTF({}, transcript:{}, start:{}, length:{})'.format(self.seqname,
                                                                        self.transcript_name,
                                                                        self.start,
                                                                        self.end-self.start+1)
        elif self.feature == 'exon':
            return 'GTF({}, exon:{}, start:{}, length:{}, id:{})'.format(self.seqname,
                                                                         self.exon_number,
                                                                         self.start,
                                                                         self.end-self.start+1,
                                                                         self.exon_id)
        elif self.feature == 'CDS':
            return 'GTF({}, CDS:{}, start:{}, length:{}, id:{})'.format(self.seqname,
                                                                         self.exon_number,
                                                                         self.start,
                                                                         self.end-self.start+1,
                                                                         self.exon_id)

    def __hash__(self):
        if self.feature == 'gene':
            return hash(self.gene_id)
        elif self.feature == 'transcript':
            return hash(self.transcript_id)
        elif self.feature == 'CDS':
            return hash('CDS' + self.exon_id)
        elif self.feature == 'exon':
            return hash(self.exon_id)

    def __eq__(self, other):
        try:
            if self.feature == 'gene':
                return self.gene_id == other.gene_id
            elif self.feature == 'transcript':
                return self.transcript_id == other.transcript_id
            elif self.feature == 'CDS':
                if other.feature == 'CDS':
                    return self.exon_id == other.exon_id
                else:
                    return False
            elif self.feature == 'exon':
                return self.exon_id == other.exon_id
            else:
                return False
        except AttributeError:
            return False