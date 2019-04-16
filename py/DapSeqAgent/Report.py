def generate_text_report(task, outfile):
    with open(outfile, 'w') as of:
        of.write('Genome\tProtein\tSample ID\tTreatment\tReplicate\tInput\tPeak ID\tSequence ID\tStart\tEnd\tSummit\tLength\tPileup\t-log(pvalue)\tEnrichment\t-log(qvalue)\tStatus\tUTRs\tGenes\tOperons\tSites\n')
        for sample in task.samples.keys():
            if 'peaks' not in task.samples[sample]:
                continue
            for peak in task.samples[sample]['peaks']:
                of.write('\t'.join([task.samples[sample]['ref_library'],
                                    task.samples[sample]['protein'],
                                    sample,
                                    task.samples[sample]['treatment'],
                                    task.samples[sample]['replicate'],
                                    task.samples[sample]['control'],
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
                                    peak.status
                                    ]))
                of.write('\t')
                if len(peak.utrs) > 0:
                    of.write(';'.join(peak.utrs))
                of.write('\t')
                if len(peak.genes) > 0:
                    of.write(';'.join(peak.genes))
                of.write('\t')
                if len(peak.operons) > 0:
                    of.write(';'.join(peak.operons))
                of.write('\t')
                if 'motifs' in task.samples[sample]:
                    sites = []
                    for motif in task.samples[sample]['motifs']:
                        for site in motif.sites:
                            if site.site_id == peak.peak_id:
                                sites.append(str(site.start) + ':' + site.sequence)
                    of.write(';'.join(sites))
                of.write('\n')
    
def generate_mapping_report(task, outfile):
    stats = ['reads processed',
            'reads aligned',
            'reads with at least one reported alignment',
            'reads with at least one reported alignment%',
            'reads that failed to align',
            'reads that failed to align%',
            'reads with alignments suppressed due to -m',
            'reads with alignments suppressed due to -m%'
            ]
    with open(outfile, 'w') as of:
        of.write('Sample\tAverage read length\tReads, total\tReads, mapped\tReads, % mapped\tReads, unaligned\tReads, % unaligned\tReads, multiple mapping\tReads, % multiple mapping\n')
        for sample_id in sorted(task.samples.keys()):
            of.write(sample_id)
            of.write('\t{0:.0f}'.format(task.samples[sample_id]['info']['average_read_length']))
            for k in stats:
                if k in task.samples[sample_id]:
                    of.write('\t' + str(task.samples[sample_id][k]))
                else:
                    of.write('\t')
            of.write('\n')
