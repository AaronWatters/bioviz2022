
"""
Parse and convert CSV challenge tiles.
Export to JSON.
"""

import json
import numpy as np

class Record:
    def __init__(self, names, values):
        name2value = {}
        for (n,v) in zip(names, values):
            name2value[n] = v
        self.names = names
        self.name2value = name2value
    def __getattr__(self, name):
        n2v = self.name2value
        if name in n2v:
            return n2v[name]
        raise AttributeError("not found: " + repr(name))
    def __repr__(self):
        inner = ["%s=%s" % (n, repr(v)) for (n,v) in sorted(self.name2value.items())]
        return "R(%s)" % ", ".join(inner)

def csv_record_list(filename):
    f = open(filename)
    headers = nextsplit(f)
    result = []
    done = False
    while not done:
        values = nextsplit(f)
        if values is None:
            done = True
        else:
            r = Record(headers, values)
            result.append(r)
    return result

class ProteinCollection:

    def __init__(self):
        self.name_to_protein = {}

    def dump(self):
        n2p = self.name_to_protein
        for name in sorted(n2p.keys()):
            p = n2p[name]
            #print(name)
            p.dump()

    def json_object(self):
        return {name: protein.json_object() for (name, protein) in self.name_to_protein.items()}

    def get_or_make_protein(self, name):
        n2p = self.name_to_protein
        if name in n2p:
            return n2p[name]
        else:
            result = Protein(name)
            n2p[name] = result
            return result

    def add_fold_record(self, record):
        name = record.UniAcc
        protein = self.get_or_make_protein(name)
        protein.add_fold_record(record)

    def add_test_record(self, record):
        ann = Annotation(record)
        protein = self.name_to_protein[ann.UniAcc]
        protein.add_annotation(ann)

class Protein:

    def __init__(self, name):
        self.name = name
        self.position_to_residue = {}
        self.gene = None
        self.species = None
        self.entry = None
        self.classifications = {}
        self.pathogenic = False
        self.locations = []

    def json_object(self):
        A = np.array(self.locations, dtype=np.float)
        m = A.min(axis=0)
        M = A.max(axis=0)
        D = M - m
        radius = np.linalg.norm(D)
        center = 0.5 * (M + m)
        result = dict(
            name=self.name, 
            gene=self.gene, 
            species=self.species, 
            entry=self.entry,
            pathogenic=self.pathogenic,
            center=center.tolist(),
            radius=float(radius),
            locations=self.locations,
            )
        result["classifications"] = sorted([[c,n] for (c,n) in self.classifications.items()])
        p2r = self.position_to_residue
        result["residues"] = {p: p2r[p].json_object() for p in sorted(p2r.keys())}
        return result

    def add_annotation(self, ann):
        pos = ann.POS
        res = self.position_to_residue[pos]
        res.add_annotation(ann)
        if self.gene is None:
            self.gene = ann.Gene
            self.species = ann.Species
            self.entry = ann.Entry
        else:
            assert self.gene == ann.Gene
            assert self.species == ann.Species
            assert self.entry == ann.Entry
        #self.classifications.add(ann.classification)
        clss = self.classifications
        cls = ann.classification
        clss[cls] = clss.get(cls, 0) + 1
        if ann.PathogenicMutation:
            self.pathogenic = True

    def dump(self):
        print()
        print("Protein: ", self.name, self.gene, self.species, self.entry)
        print("   ", sorted(self.classifications.items()))
        p2r = self.position_to_residue
        for p in sorted(p2r.keys()):
            r = p2r[p]
            r.dump()

    def add_fold_record(self, record):
        res = Residue(record)
        pos = res.POS
        p2r = self.position_to_residue
        assert pos not in p2r, "Duplicate residue position: " + repr((pos, list(p2r.keys())))
        p2r[pos] = res
        self.locations.append(res.location())

class Residue:

    def __init__(self, record):
        self.POS = int(record.POS)
        self.RES = str(record.RES)
        self.quality = float(record.quality)
        C = Atom("c", "C", record)
        assert C.valid
        self.C = C
        N = Atom("n", "N", record)
        assert N.valid
        self.N = N
        CA = Atom("ca", "C", record)
        assert CA.valid
        self.CA = CA
        CB = Atom("cb", "C", record)
        self.CB = None
        if CB.valid:
            self.CB = CB
        self.annotations = []
        self.classifications = {}
        self.PathogenicMutation = False

    def add_classification(self, cls):
        clss = self.classifications
        clss[cls] = clss.get(cls, 0) + 1

    def location(self):
        return self.C.location()

    def json_object(self):
        D = dict(
            POS = self.POS,
            RES = self.RES,
            quality = self.quality,
            C = self.C.json_object(),
            N = self.N.json_object(),
            CA = self.CA.json_object(),
            CB = None,
            annotations = [ann.json_object() for ann in self.annotations],
            pathogenic = self.PathogenicMutation,
        )
        if self.CB is not None:
            D["CB"] = self.CB.json_object()
        return D

    def add_annotation(self, ann):
        assert self.RES == ann.RES
        assert self.POS == ann.POS
        self.annotations.append(ann)
        if ann.PathogenicMutation:
            self.PathogenicMutation = True
        self.add_classification(ann.classification)

    def dump(self, indent="   "):
        print(indent, "RESIDUE", self.POS, self.RES, self.quality, sorted(self.classifications.items()))
        if self.PathogenicMutation:
            print(indent, "  PathogenicMutation")
        print(indent, "       ", self.C, self.CA, self.N, self.CB)
        for ann in self.annotations:
            print(indent, "    Annotation:", ann)

class Atom:

    coord_prefixes = "x_coord_ y_coord_ z_coord_".split()

    def __init__(self, identifier, element, record):
        self.valid = True
        self.id = identifier
        self.element = element
        for prefix in self.coord_prefixes:
            attr = prefix + identifier
            try:
                svalue = getattr(record, attr)
            except AttributeError:
                self.valid = False
            else:
                coord = prefix[0]
                try:
                    val = float(svalue)
                except ValueError:
                    self.valid = False
                else:
                    setattr(self, coord, val)

    def json_object(self):
        D = dict(
            id=self.id,
            element=self.element,
            pos=self.location(),
        )
        return D

    def location(self):
        return [self.x, self.y, self.z]

    def __repr__(self):
        return self.id + ":" + self.element + repr([self.x, self.y, self.z])

class Annotation:

    def __init__(self, record):
        self.UniAcc = record.UniAcc
        self.RES = record.RES
        self.POS = int(record.POS)
        self.MOD = record.MOD
        self.Entry = record.Entry
        self.Gene = record.Gene
        self.Species = record.Species
        self.classification = record.classification
        pm = record.PathogenicMutation
        p = True
        if pm.upper() == "FALSE":
            p = False
        #print("Annotation pathogenic", repr(pm), p)
        self.PathogenicMutation = p

    def __repr__(self):
        pmarker = ""
        if self.PathogenicMutation:
            pmarker = " PATHOGENIC"
        return self.MOD + ":" + self.classification + pmarker

    def json_object(self):
        return dict(MOD=self.MOD, classification=self.classification)

def nextsplit(f):
    line = f.readline()
    sline = line.strip()
    if not sline:
        return None
    ssline = sline.split(",")
    return ssline

fold_data = "BioVis-challenge-alphafold-data.csv"
test_data = "BioVis-challenge-test-data.csv"
json_out = "challenge.json"

def generate_json(fold_file_name=fold_data, test_file_name=test_data, out_fn=json_out, verbose=True):
    fold_records = csv_record_list(fold_file_name)
    #if verbose:
    #    print("========= Fold records:")
    #    for r in fold_records:
    #        print(r)
    test_records = csv_record_list(test_file_name)
    #if verbose:
    #    print("========== Test records:")
    #    for r in test_records:
    #        print(r)
    col = ProteinCollection()
    for r in fold_records:
        col.add_fold_record(r)
    for r in test_records:
        col.add_test_record(r)
    if verbose:
        print()
        print("======== COLLECTION")
        col.dump()
    json_ob = col.json_object()
    f = open(out_fn, "w")
    json.dump(json_ob, f, indent=2)
    f.close()
    print("Dumped JSON to", repr(out_fn))

if __name__ == "__main__":
    generate_json()
