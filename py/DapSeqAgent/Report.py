def generate_text_report(task, outfile):
    with open(outfile, 'w') as of:
        of.write('Genome\tProtein\tSample ID\tTreatment\tReplicate\tPeak ID\tSequence ID\tStart\tEnd\tSummit\tLength\tPileup\t-log(pvalue)\tEnrichment\t-log(qvalue)\tSummit in UTR\tUTRs\tUTR functions\tGenes\tGene functions\tOperon\n')
        for sample in task.samples.keys():
            if 'peaks' not in task.samples[sample]:
                continue
            for peak in task.samples[sample]['peaks']:
                of.write('\t'.join([task.samples[sample]['ref_library'],
                                    task.samples[sample]['protein'],
                                    sample,
                                    task.samples[sample]['treatment'],
                                    task.samples[sample]['replicate'],
                                    peak.peak_id,
                                    peak.sequence_id,
                                    str(peak.start),
                                    str(peak.end),
                                    str(peak.abs_summit),
                                    str(peak.length),
                                    str(peak.pileup),
                                    str(peak.log_pvalue),
                                    str(peak.fold_enrichment),
                                    str(peak.log_qvalue),
                                    str(peak.noncoding)
                                    ]))
                of.write('\t')
                if len(peak.utrs) > 0:
                    of.write(';'.join(peak.utrs))
                    of.write('\t')
                    of.write(peak.utr_functions)
                else:
                    of.write('\t')
                of.write('\t')
                if len(peak.genes) > 0:
                    of.write(';'.join(peak.genes))
                    of.write('\t')
                    of.write(peak.gene_functions)
                else:
                    of.write('\t')
                of.write('\t')
                if len(peak.operons) > 0:
                    of.write(';'.join(peak.operons))
                of.write('\n')
    
