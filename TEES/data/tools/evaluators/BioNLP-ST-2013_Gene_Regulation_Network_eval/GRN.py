#!/usr/bin/python


# Copyright (c) 2013, Institut National de la Recherche Agronomique (INRA)
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

#     Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#     Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
#     Neither the names of the Institut National de la Recherche Agronomique (INRA) and BioNLP-ST 2013 nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



from sys import exit
from bionlpst import *
from os import mkdir
from os.path import exists
import re
from StringIO import StringIO

GRN_SCHEMA = '''[DEFAULT]
interaction_agent = Gene, GeneFamily, mRNA, Protein, ProteinComplex, PolymeraseComplex, ProteinFamily, Action_Target
interaction_target = Gene, Operon, GeneFamily, mRNA, Protein, ProteinComplex, PolymeraseComplex, ProteinFamily, Action_Target, Transcription_by, Transcription_from, Interaction.Transcription, Interaction.Activation

[Action]
kind = entity

[Protein]
kind = entity

[Gene]
kind = entity

[PolymeraseComplex]
kind = entity

[Site]
kind = entity

[GeneFamily]
kind = entity

[Promoter]
kind = entity

[ProteinFamily]
kind = entity

[Operon]
kind = entity

[Regulon]
kind = entity

[ProteinComplex]
kind = entity

[mRNA]
kind = entity

[Gene_Identifier]
kind = normalization
target = Gene, GeneFamily, Operon, Protein, ProteinComplex, PolymeraseComplex, ProteinFamily, mRNA
referent = .*

[Action_Target]
kind = event
trigger = Action
arg.Target = Gene, Operon, GeneFamily, mRNA, Protein, ProteinComplex, ProteinFamily, PolymeraseComplex, Promoter, Regulon

[Transcription_by]
kind = event
trigger = Action
arg.Agent = PolymeraseComplex, Protein, ProteinComplex, ProteinFamily, Gene, GeneFamily

[Transcription_from]
kind = event
trigger = Action
arg.Site = Promoter, Site

[Site_of]
kind = relation
arg.Entity = Gene, Operon, Promoter
arg.Site = Site, Promoter

[Member_of_Regulon]
kind = relation
arg.Regulon = Regulon
arg.Member = Gene, Operon, Protein, ProteinComplex

[Master_of_Regulon]
kind = relation
arg.Regulon = Regulon
arg.Master = Gene, Protein

[Promoter_of]
kind = relation
arg.Promoter = Promoter
arg.Gene = Gene, Operon

[Master_of_Promoter]
kind = relation
arg.Promoter = Promoter
arg.Protein = Protein, ProteinComplex, PolymeraseComplex, Gene

[Bind_to]
kind = relation
arg.DNA = Gene, Site
arg.Protein = Protein, PolymeraseComplex, ProteinComplex

[Interaction.Regulation]
kind = relation
arg.Agent = %(interaction_agent)s
arg.Target = %(interaction_target)s

[Interaction.Activation]
kind = relation
arg.Agent = %(interaction_agent)s
arg.Target = %(interaction_target)s

[Interaction.Inhibition]
kind = relation
arg.Agent = %(interaction_agent)s
arg.Target = %(interaction_target)s

[Interaction.Requirement]
kind = relation
arg.Agent = %(interaction_agent)s
arg.Target = %(interaction_target)s

[Interaction.Transcription]
kind = relation
arg.Agent = %(interaction_agent)s
arg.Target = %(interaction_target)s

[Interaction.Binding]
kind = relation
arg.Agent = %(interaction_agent)s
arg.Target = %(interaction_target)s

[Negation]
kind = modification
target = Interaction.Regulation, Interaction.Activation, Interaction.Inhibition, Interaction.Requirement, Interaction.Binding, Interaction.Transcription
'''

def gene_id(annotation):
    for gid in annotation.get_norms('Gene_Identifier'):
        return gid
    return None

def negated(annotation):
    for neg in annotation.get_mods('Negation'):
        return True
    return False

def nodes(annotation, is_agent):
    if isinstance(annotation, TextBound):
        gid = gene_id(annotation)
        if gid is not None:
            yield annotation
        for a in annotation.doc.annotations.itervalues():
            for n in _bounce_nodes(a, annotation, is_agent):
                yield n
    if negated(annotation):
        raise Exception(annotation.message('negated interaction argument'))
    if isinstance(annotation, AnnotationWithArgs):
        for a in annotation.args.itervalues():
            for n in nodes(a, is_agent):
                yield n
        if isinstance(annotation, Event):
            for n in nodes(annotation.trigger, is_agent):
                yield n

def _bounce_nodes(annotation, arg, is_agent):
    if annotation.type == 'Promoter_of' and annotation.args['Promoter'] == arg and not is_agent:
        for n in nodes(annotation.args['Gene'], is_agent):
            yield n
    if annotation.type == 'Master_of_Promoter' and annotation.args['Promoter'] == arg and is_agent:
        for n in nodes(annotation.args['Protein'], is_agent):
            yield n

def arcs(annotation):
    if annotation.type.startswith('Interaction.') and not negated(annotation):
        agents = tuple(nodes(annotation.args['Agent'], True))
        targets = tuple(nodes(annotation.args['Target'], False))
        if len(agents) == 0:
            yield annotation.message('agent could not be associated to a gene (target: ' + ', '.join(str(t) for t in targets) + ')\n')
        if len(targets) == 0:
            yield annotation.message('target could not be associated to a gene (agent: ' + ', '.join(str(t) for t in agents) + ')\n')
        next = max(int(r.id[1:]) for r in annotation.doc.annotations.itervalues() if isinstance(r, Relation))
        for a in agents:
            for t in targets:
                next += 1
                rel = Relation(annotation.source, -1, annotation.doc, annotation.visibility, 'R' + str(next), annotation.type, {'Agent': a.id, 'Target': t.id})
                rel.original = annotation
                rel.resolve_ids()

INTERACTION_PRECEDENCE = {
    'Interaction.Regulation': (),
    'Interaction.Binding': set(('Interaction.Regulation',)),
    'Interaction.Transcription': set(('Interaction.Binding', 'Interaction.Regulation')),
    'Interaction.Other': set(('Interaction.Regulation',)),
    'Interaction.Inhibition': set(('Interaction.Regulation',)),
    'Interaction.Activation': set(('Interaction.Regulation',)),
    'Interaction.Requirement': set(('Interaction.Regulation','Interaction.Activation')),
    }


def is_specialization(a, l):
    for b in l:
        try:
            if a in INTERACTION_PRECEDENCE[b]:
                return b
        except KeyError:
            pass
    return False



def resolve_arc(arcs):
    if 'Regulation' in arcs:
        return 'Regulation'
    return get_mechanism(arcs), get_level(arcs)

MECHANISMS = set(('Interction.Transcription', 'Interaction.Binding'))
def get_mechanism(arcs):
    for m in MECHANISMS:
        if m in arcs:
            return m
    return None

LEVELS = set(('Interaction.Inhibition', 'Interaction.Activation', 'Interaction.Requirement'))
def get_level(arcs):
    for m in LEVELS:
        if m in arcs:
            return m
    return None



class Evaluation:
    def __init__(self):
        self.ref = []
        self.pred = []
        self.ins = []
        self.dels = []
        self.sub = []
        self.match = []

    def _evaluate_pair(self, pair, ref, pred):
        for r in ref:
            self.ref.append((pair, r))
        for p in pred:
            self.pred.append((pair, p))
        for a in set(ref) & set(pred):
            self.match.append((pair, a))
        ins0 = set(pred) - set(ref)
        del0 = set(ref) - set(pred)
        for i in tuple(ins0):
            d = is_specialization(i, del0)
            if d:
                self.sub.append((pair, (d, i)))
                ins0.remove(i)
                del0.remove(d)
        for d in tuple(del0):
            i = is_specialization(d, ins0)
            if i:
                self.sub.append((pair, (d, i)))
                ins0.remove(i)
                del0.remove(d)
        for i in tuple(ins0):
            for d in tuple(del0):
                self.sub.append((pair, (d, i)))
                ins0.remove(i)
                del0.remove(d)
                break
        for i in ins0:
            self.ins.append((pair, i))
        for d in del0:
            self.dels.append((pair, d))


class BioNLP_ST_GRN(BioNLP_ST):
    def __init__(self):
        BioNLP_ST.__init__(self, StringIO(GRN_SCHEMA))
        self.add_option('--interactions-dir', action='store', type='string', dest='interactions_dir', help='directory where to write resolved interactions in BioNLP-ST format')
        self.add_option('--sif', action='store', type='string', dest='sif', help='name of the file where to write the reference regulation network in Cytoscape\'s SIF format')
        self.add_option('--pred-dir', action='store', type='string', dest='pred_dir', help='name of the directory containing prediction .a2 files')
        self.add_option('--pred-sif', action='store', type='string', dest='pred_sif', help='name of the file containing the predicted regulation network in Cytoscape\'s SIF format (Agent, Interaction, Target)')
        self.add_option('--verbose', action='store_true', dest='verbose', help='print detailed error analysis')

    def _compute_arcs(self, corpus):
        nerr = 0
        print 'Resolving interactions:'
        for a in tuple(corpus.iterannotations()):
            for msg in arcs(a):
                stderr.write(msg)
                nerr += 1
        if nerr == 0:
            print '    ok'
        else:
            print '   ', str(nerr), 'errors'

    def _write_resolved_interactions(self, corpus):
        if self.options.interactions_dir is None:
            return
        print 'Writing resolved interactions in BioNLP-ST format:', self.options.interactions_dir
        if not exists(self.options.interactions_dir):
            mkdir(self.options.interactions_dir)
        for doc in corpus.iterdocs():
            fn = self.options.interactions_dir + '/' + doc.id + '.a2'
            f = open(fn, 'w')
            for a in doc.iterannotations():
                if a.lineno == -1:
                    f.write('%s\t%s Target:%s Agent:%s\n' % (a.id, a.type, a.args['Target'].id, a.args['Agent'].id))
            f.close()

    def _get_arcs_from_corpus(self, corpus):
        result = {}
        for a in corpus.iterannotations():
            if a.lineno == -1:
                pair = gene_id(a.args['Agent']).referent, gene_id(a.args['Target']).referent
                if pair in result:
                    arc = result[pair]
                else:
                    arc = {}
                    result[pair] = arc
                if a.type in arc:
                    v = arc[a.type]
                else:
                    v = []
                    arc[a.type] = v
                v.append(a.original)
        return self._remove_all_redundant(result)

    def _get_arcs_from_sif(self, source):
        result = {}
        lineno = 0
        file = open(source)
        for line in file:
            lineno += 1
            line = line.strip()
            if line == '':
                continue
            agent, inter, target = line.split('\t')
            pair = agent, target
            if pair in result:
                arc = result[pair]
            else:
                arc = {}
                result[pair] = arc
            if inter in arc:
                v = arc[a.type]
            else:
                v = []
                arc[inter] = v
            v.append(Sourced(source, lineno, None, A2))
        return self._remove_all_redundant(result)

    def _remove_redundant(self, arc):
        all_types = arc.keys()
        return dict((type, v) for (type, v) in arc.iteritems() if not is_specialization(type, all_types))

    def _remove_all_redundant(self, arcs):
        return dict((pair, self._remove_redundant(arc)) for (pair, arc) in arcs.iteritems())
   
    def _write_sif(self, corpus, sif):
        if sif is None:
            return
        all_arcs = self._get_arcs_from_corpus(corpus)
        print 'Writing resolved interactions in SIF (cytoscape) format:', sif
        f = open(sif, 'w')
        for (agent, target), arcs in all_arcs.iteritems():
            for arc in arcs:
                if not is_specialization(arc, arcs):
                    f.write('%s\t%s\t%s\n' % (agent, arc, target))
        f.close()

    def run(self):
        BioNLP_ST.run(self)
        self._compute_arcs(self.corpus)
        self._write_resolved_interactions(self.corpus)
        self._write_sif(self.corpus, self.options.sif)
        if self.options.pred_dir or self.options.pred_sif:
            if self.options.pred_dir:
                self.prediction = self._load_corpus(self.options.pred_dir)
                self._compute_arcs(self.prediction)
            self._evaluate()
            self._write_evaluation()

    def _evaluate(self):
        self.evaluate = Evaluation()
        ref = self._get_arcs_from_corpus(self.corpus)
        if self.options.pred_sif is not None:
            pred = self._get_arcs_from_sif(self.options.pred_sif)
        else:
            pred = self._get_arcs_from_corpus(self.prediction)
        for pair, ra in ref.iteritems():
            if pair in pred:
                self.evaluate._evaluate_pair(pair, ra, pred[pair])
            else:
                self.evaluate._evaluate_pair(pair, ra, {})
        for pair, pa in pred.iteritems():
            if pair not in ref:
                self.evaluate._evaluate_pair(pair, {}, pa)

    def _write_evaluation(self):
        print
        print 'Evaluation:'
        print
        nsub = len(self.evaluate.sub)
        ndel = len(self.evaluate.dels)
        nins = len(self.evaluate.ins)
        nok = len(self.evaluate.match)
        nref = len(self.evaluate.ref)
        npred = len(self.evaluate.pred)
        print 'Substitutions:', nsub
        print 'Deletions:', ndel
        print 'Insertions:', nins
        print 'Matches:', nok
        print 'Reference:', nref
        print 'Predictions:', npred
        print
        print 'Slot Error Rate:', (float(nsub + ndel + nins) / float(nref))
        r = float(nok) / float(nref)
        print
        print 'Recall:', r
        p = float(nok) / float(npred)
        print 'Precision:', p
        print 'F-score:', ((2 * r * p) / (r + p))
        print
        print 'Relaxed Slot Error Rate:', (float(ndel + nins) / float(nref))
        r = float(nok + nsub) / float(nref)
        print
        print 'Relaxed Recall:', r
        p = float(nok + nsub) / float(npred)
        print 'Relaxed Precision:', p
        if r + p == 0:
            f = 0
        else:
            f = ((2 * r * p) / (r + p))
        print 'Relaxed F-score:', f
        print
        if (self.options.verbose):
            ref = self._get_arcs_from_corpus(self.corpus)
            if self.options.pred_sif is not None:
                pred = self._get_arcs_from_sif(self.options.pred_sif)
            else:
                pred = self._get_arcs_from_corpus(self.prediction)
#            pred = self._get_arcs(self.prediction)
            print
            print 'Details:'
            print
            print 'Substitutions:'
            for (a, t), (r, p) in self.evaluate.sub:
                print '    %s/%s: %s -> %s' % (a, t, r, p)
                print '        Reference sources:'
                for ann in ref[a, t][r]:
                    print '           ', ann.message(ann.id)
                print '        Prediction sources:'
                for ann in pred[a, t][p]:
                    print '           ', ann.message(ann.id)
            print
            print 'Deletions:'
            for (a, t), r in self.evaluate.dels:
                print '    %s/%s: %s' % (a, t, r)
                print '        Reference sources:'
                for ann in ref[a, t][r]:
                    print '           ', ann.message(ann.id)
            print
            print 'Insertions:'
            for (a, t), p in self.evaluate.ins:
                print '    %s/%s: %s' % (a, t, r)
                print '        Prediction sources:'
                for ann in pred[a, t][p]:
                    print '           ', ann.message(ann.id)
            print
        
    
if __name__ == '__main__':
    BioNLP_ST_GRN().run()
