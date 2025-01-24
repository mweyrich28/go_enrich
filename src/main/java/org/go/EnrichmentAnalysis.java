package org.go;

import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.go.util.FileUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class EnrichmentAnalysis
{
    private static final Logger logger = LoggerFactory.getLogger(EnrichmentAnalysis.class);

    private DAG dag;
    private HashMap<String, Gene> enrichedGeneMap;
    private String overlapOut;
    private String out;
    private int min;
    private int max;
    private HashSet<String> enrichedSigs = new HashSet<>();
    private HashSet<String> enrichedAll;

    public EnrichmentAnalysis(DAG dag, HashMap<String, Gene> enrichedGeneMap, int min, int max, String overlapOut, String out)
    {
        this.dag = dag;
        this.enrichedGeneMap = enrichedGeneMap;
        this.overlapOut = overlapOut;
        this.out = out;
        this.min = min;
        this.max = max;

        this.enrichedAll = new HashSet<>(enrichedGeneMap.keySet());
        for (Gene gene : enrichedGeneMap.values())
        {
            if (gene.isSig())
            {
                this.enrichedSigs.add(gene.getSymbol());
            }
        }
    }

    public void analyze()
    {
        logger.info("Starting gene set enrichment analysis...");
        double start = System.currentTimeMillis();
        StringBuilder sb = new StringBuilder();
        sb.append("term\tname\tsize\tis_true\tnoverlap\thg_pval\thg_fdr\tfej_pval\tfej_fdr\tks_stat\tks_pval\tks_fdr\tshortest_path_to_a_true\n");

        for (GOEntry go : dag.getNodeMap().values())
        {
            if (!(go.getGeneSymbols().size() >= min && go.getGeneSymbols().size() <= max))
            {
                continue;
            }

            // id
            sb.append(go.getId()).append("\t");
            // annot
            sb.append(go.getAnnot()).append("\t");

            // size
            HashSet<String> intersection = intersect(go, enrichedAll);
            sb.append(intersection.size()).append("\t");

            // is_true
            sb.append(go.isTrue()).append("\t");

            // noverlap
            intersection = intersect(go, enrichedSigs);
            sb.append(intersection.size()).append("\t");

            System.out.println(sb.toString());
            sb.setLength(0);

        }
        logger.info(String.format("Time needed for GSEA: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
    }

    public HashSet<String> intersect(GOEntry go, HashSet<String> set)
    {
        HashSet<String> intersection = new HashSet<>(go.getGeneSymbols());
        intersection.retainAll(set);
        return intersection;
    }
}
