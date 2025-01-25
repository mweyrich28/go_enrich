package org.go;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
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
    private HashSet<String> enrichedNotSigs = new HashSet<>();

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
            else
            {
                enrichedNotSigs.add(gene.getSymbol());
            }
        }
    }

    public void analyze() throws IOException
    {
        logger.info("Starting gene set enrichment analysis...");
        BufferedWriter buff = new BufferedWriter(new FileWriter(out));
        double start = System.currentTimeMillis();
        //        sb.append("term\tname\tsize\tis_true\tnoverlap\thg_pval\thg_fdr\tfej_pval\tfej_fdr\tks_stat\tks_pval\tks_fdr\tshortest_path_to_a_true");
        ArrayList<Double> hgPvalsList = new ArrayList<>();
        ArrayList<Double> fejPvalsList = new ArrayList<>();
        ArrayList<AnalysisEntry> rows = new ArrayList<>();

        for (GOEntry go : dag.getNodeMap().values())
        {
            if (!(go.getGeneSymbols().size() >= min && go.getGeneSymbols().size() <= max))
            {
                continue;
            }

            // -------------------- id -------------------------
            AnalysisEntry entry = new AnalysisEntry(go.getId());
            rows.add(entry);

            // --------- annot ----------
            entry.setName(go.getAnnot());

            // -------------------------------- size -----------------------------------------
            HashSet<String> overlappingEnriched = intersect(go.getGeneSymbols(), enrichedAll);
            entry.setSize(overlappingEnriched.size());

            // ------- is_true -----------
            entry.setIs_true(go.isTrue());

            // ----------------------- noverlap -----------------------------------------------
            HashSet<String> overlappingSigsInGo = intersect(go.getGeneSymbols(), enrichedSigs);
            entry.setNoverlap(overlappingSigsInGo.size());

            // --------------------------- hg_pval ---------------------------------------------------
            HashSet<String> overlappingTotal = intersect(dag.getRoot().getGeneSymbols(), enrichedAll);
            HashSet<String> overlappingSigs = intersect(dag.getRoot().getGeneSymbols(), enrichedSigs);
            double hgpval = calcHgPval(overlappingTotal, overlappingEnriched, overlappingSigs, overlappingSigsInGo);
            entry.setHg_pval(hgpval);
            hgPvalsList.add(hgpval);

            // --------- fej_pval -------------------
            double fejpval = jackknifeHgPval(overlappingTotal, overlappingEnriched, overlappingSigs, overlappingSigsInGo);
            entry.setFej_pval(fejpval);
            fejPvalsList.add(fejpval);

            // ----------- ks_stat ------------------
            // ----------- ks_pval ------------------
            // ----------- ks_fdr -------------------
            // ----------- shortest_path_to_a_true ------------------


            //            break;
        }
        // -------------------- hg_fdr -----------------------------
        ArrayList<Double> adjHgPvals = calculateBH_FDR(hgPvalsList);
        ArrayList<Double> adjFejPvals = calculateBH_FDR(fejPvalsList);
        int i = 0;
        for (AnalysisEntry entry: rows)
        {
            entry.setHg_fdr(adjHgPvals.get(i));
            entry.setFej_fdr(adjFejPvals.get(i));
            i++;
            System.out.println(entry);
            if (i == 4) break;
        }


        // ----------- fej_fdr ------------------

        buff.flush();
        logger.info(String.format("Time needed for Analysis: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
    }

    public HashSet<String> intersect(HashSet<String> A, HashSet<String> B)
    {
        HashSet<String> intersection = new HashSet<>(A);
        intersection.retainAll(B);
        return intersection;
    }

    public double calcHgPval(HashSet<String> overlappingTotal,
                             HashSet<String> overlappingEnriched,
                             HashSet<String> overlappingSigs,
                             HashSet<String> overlappingSigsInGO)
    {
        int N = overlappingTotal.size();
        int n = overlappingEnriched.size();
        int K = overlappingSigs.size();
        int k = overlappingSigsInGO.size();

        HypergeometricDistribution hg = new HypergeometricDistribution(N, K, n);
        return hg.upperCumulativeProbability(k);
    }


    public ArrayList<Double> calculateBH_FDR(ArrayList<Double> pValues) {
        int m = pValues.size();
        ArrayList<Integer> ranks = new ArrayList<>();
        ArrayList<Double> sortedPValues = new ArrayList<>(pValues);
        ArrayList<Double> bhAdjustedPValues = new ArrayList<>(Collections.nCopies(m, 0.0));

        HashMap<Double, Integer> pValueToIndex = new HashMap<>();
        for (int i = 0; i < m; i++) {
            pValueToIndex.put(pValues.get(i), i);
        }

        sortedPValues.sort(Double::compareTo);

        for (int i = 0; i < m; i++) {
            int rank = i + 1; // Rank starts at 1
            double rawPValue = sortedPValues.get(i);
            double adjustedPValue = rawPValue * m / rank;
            bhAdjustedPValues.set(i, Math.min(adjustedPValue, 1.0));
        }

        for (int i = m - 2; i >= 0; i--) {
            if (bhAdjustedPValues.get(i) > bhAdjustedPValues.get(i + 1)) {
                bhAdjustedPValues.set(i, bhAdjustedPValues.get(i + 1));
            }
        }

        ArrayList<Double> finalAdjustedPValues = new ArrayList<>(Collections.nCopies(m, 0.0));
        for (int i = 0; i < m; i++) {
            int originalIndex = pValueToIndex.get(sortedPValues.get(i));
            finalAdjustedPValues.set(originalIndex, bhAdjustedPValues.get(i));
        }

        return finalAdjustedPValues;
    }

    public double jackknifeHgPval(HashSet<String> overlappingTotal,
                                             HashSet<String> overlappingEnriched,
                                             HashSet<String> overlappingSigs,
                                             HashSet<String> overlappingSigsInGO) {
        ArrayList<Double> jackknifePvals = new ArrayList<>();

        ArrayList<String> sigsInGOList = new ArrayList<>(overlappingSigsInGO);

        for (int i = 0; i < sigsInGOList.size(); i++) {
            HashSet<String> jackknifeSigsInGO = new HashSet<>(overlappingSigsInGO);
            jackknifeSigsInGO.remove(sigsInGOList.get(i));

            int N = overlappingTotal.size();
            int n = overlappingEnriched.size();
            int K = overlappingSigs.size();
            int k = jackknifeSigsInGO.size();

            HypergeometricDistribution hg = new HypergeometricDistribution(N, K, n);
            double pval = hg.upperCumulativeProbability(k);

            jackknifePvals.add(pval);
        }

        double meanPval = jackknifePvals.stream()
                .mapToDouble(Double::doubleValue)
                .average()
                .orElse(0.0); // Handle empty list case with 0.0 as default
        return meanPval;
    }
}
