package org.go;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.*;
import java.util.stream.Collectors;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class EnrichmentAnalysis
{

    private static final Logger logger = LoggerFactory.getLogger(EnrichmentAnalysis.class);
    private DAG dag;
    private HashMap<String, Gene> enrichedGeneMap;
    private KolmogorovSmirnovTest ks = new KolmogorovSmirnovTest();
    private String out;
    private int min;
    private int max;
    private HashSet<String> enrichedSigs = new HashSet<>();
    private HashSet<String> enrichedAll;
    private HashSet<String> enrichedNotSigs = new HashSet<>();

    public EnrichmentAnalysis(DAG dag, HashMap<String, Gene> enrichedGeneMap, int min, int max, String out)
    {

        this.dag = dag;
        this.enrichedGeneMap = enrichedGeneMap;
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
                this.enrichedNotSigs.add(gene.getSymbol());
            }
        }
    }

    public void goFeatures(String pathToOverlapOut) throws IOException
    {
        logger.info("Collecting GO features...");
        BufferedWriter buff = new BufferedWriter(new FileWriter(pathToOverlapOut));
        buff.write("term1\tterm2\tis_relative\tpath_length\tnum_overlapping\tmax_ov_percent");
        double start = System.currentTimeMillis();

        ArrayList<GOEntry> goEntries = new ArrayList<>(dag.getNodeMap().values());
        for (int i = 0; i < goEntries.size(); i++)
        {
            GOEntry term1 = goEntries.get(i);
            if (!(term1.getGeneSymbols().size() >= min && term1.getGeneSymbols().size() <= max))
            {
                continue;
            }
            for (int j = i + 1; j < goEntries.size(); j++)
            {
                GOEntry term2 = goEntries.get(j);
                if (!(term2.getGeneSymbols().size() >= min && term2.getGeneSymbols().size() <= max))
                {
                    continue;
                }

                int numOverlapping = smartIntersect(term1.getGeneSymbols(), term2.getGeneSymbols());

                if (numOverlapping == 0)
                {
                    continue;
                }

                boolean isRealtive = false;
                if (term1.getReachableGoIDs().contains(term2.getId()) || term2.getReachableGoIDs().contains(term1.getId()))
                {
                    isRealtive = true;
                }

                // moved logic to extra class, too crowded in DAG and EA
                PathComputer smartPath = new PathComputer(term1, term2);
                int shortestPath = smartPath.calculateShortestPath();
                double maxOvPct = 100.0 * (double) numOverlapping / Math.min(term1.getGeneSymbols().size(), term2.getGeneSymbols().size());

                buff.write("\n" + term1.getId() + "\t" + term2.getId() + "\t" + isRealtive + "\t" + shortestPath + "\t" + numOverlapping + "\t" + maxOvPct);
            }
        }
        buff.flush();
        logger.info(String.format("Time needed for Collecting GO features: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
    }

    public void goFeaturesParallelized(String pathToOverlapOut) throws IOException
    {
        logger.info("Collecting GO features...");
        double start = System.currentTimeMillis();

        List<GOEntry> goEntries = new ArrayList<>(dag.getNodeMap().values());
        int size = goEntries.size();

        ConcurrentLinkedQueue<String> results = new ConcurrentLinkedQueue<>();
        results.add("term1\tterm2\tis_relative\tpath_length\tnum_overlapping\tmax_ov_percent");

        IntStream.range(0, size).parallel().forEach(i -> {
            GOEntry term1 = goEntries.get(i);
            if (!(term1.getGeneSymbols().size() >= min && term1.getGeneSymbols().size() <= max))
            {
                return;
            }

            for (int j = i + 1; j < size; j++)
            {
                GOEntry term2 = goEntries.get(j);
                if (!(term2.getGeneSymbols().size() >= min && term2.getGeneSymbols().size() <= max))
                {
                    continue;
                }

                int numOverlapping = smartIntersect(term1.getGeneSymbols(), term2.getGeneSymbols());
                if (numOverlapping == 0)
                {
                    continue;
                }

                boolean isRelative = term1.getReachableGoIDs().contains(term2.getId()) || term2.getReachableGoIDs().contains(term1.getId());

                PathComputer smartPath = new PathComputer(term1, term2);
                int shortestPath = smartPath.calculateShortestPath();

                double maxOvPct = 100.0 * (double) numOverlapping / Math.min(term1.getGeneSymbols().size(), term2.getGeneSymbols().size());

                results.add(term1.getId() + "\t" + term2.getId() + "\t" + isRelative + "\t" + shortestPath + "\t" + numOverlapping + "\t" + maxOvPct);
            }
        });

        BufferedWriter buff = new BufferedWriter(new FileWriter(pathToOverlapOut));
        for (String line : results)
        {
            buff.write(line);
            buff.write("\n");
        }
        buff.flush();
        logger.info(String.format("Time needed for Collecting GO features: %s seconds", (System.currentTimeMillis() - start) / 1000.0));

    }

    public void analyze() throws IOException
    {
        logger.info("Starting gene set enrichment analysis...");
        BufferedWriter buff = new BufferedWriter(new FileWriter(out));
        double start = System.currentTimeMillis();
        buff.write("term\tname\tsize\tis_true\tnoverlap\thg_pval\thg_fdr\tfej_pval\tfej_fdr\tks_stat\tks_pval\tks_fdr\tshortest_path_to_a_true");
        ArrayList<Double> hgPvalsList = new ArrayList<>();
        ArrayList<Double> fejPvalsList = new ArrayList<>();
        ArrayList<Double> ksPvalsList = new ArrayList<>();
        ArrayList<AnalysisEntry> rows = new ArrayList<>();

        for (GOEntry go : dag.getNodeMap().values())
        {
            if (!(go.getGeneSymbols().size() >= min && go.getGeneSymbols().size() <= max))
            {
                continue;
            }

            AnalysisEntry entry = new AnalysisEntry(go.getId());
            rows.add(entry);


            entry.setName(go.getAnnot());

            int n = smartIntersect(go.getGeneSymbols(), enrichedAll);
            entry.setSize(n);

            entry.setIs_true(go.isTrue());

            int k = smartIntersect(go.getGeneSymbols(), enrichedSigs);
            entry.setNoverlap(k);

            int N = smartIntersect(dag.getRoot().getGeneSymbols(), enrichedAll);
            int K = smartIntersect(dag.getRoot().getGeneSymbols(), enrichedSigs);
            double hgpval = calcHgPval(N, n, K, k);
            entry.setHg_pval(hgpval);
            hgPvalsList.add(hgpval);

            double fejpval = calFejPval(N, n, K, k);
            entry.setFej_pval(fejpval);
            fejPvalsList.add(fejpval);

            HashSet<String> bg = new HashSet<>(dag.getRoot().getGeneSymbols());
            HashSet<String> overlappingEnriched = intersect(go.getGeneSymbols(), enrichedAll);
            bg.removeAll(overlappingEnriched);

            double[] ksStats = calcKS(overlappingEnriched, bg);
            entry.setKs_stat(ksStats[0]);
            entry.setKs_pval(ksStats[1]);
            ksPvalsList.add(ksStats[1]);

            if (go.isTrue() || dag.getTrueGoEntries().isEmpty())
            {
                entry.setShortest_path_to_a_true("");
            }
            else
            {
                String path = go.getShortestPathToTrue(dag.getTrueGoIds(), dag);
                entry.setShortest_path_to_a_true(path);
            }

        }

        ArrayList<Double> adjHgPvals = calculateBH_FDR(hgPvalsList);
        ArrayList<Double> adjFejPvals = calculateBH_FDR(fejPvalsList);
        ArrayList<Double> adjKsPvals = calculateBH_FDR(ksPvalsList);
        int i = 0;
        for (AnalysisEntry entry : rows)
        {

            entry.setHg_fdr(adjHgPvals.get(i));
            entry.setFej_fdr(adjFejPvals.get(i));
            entry.setKs_fdr(adjKsPvals.get(i));
            i++;
            buff.write("\n");
            buff.write(entry.toString());
        }


        buff.flush();
        logger.info(String.format("Time needed for Analysis: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
    }

    public void analyzeParallelized() throws IOException
    {
        logger.info("Starting gene set enrichment analysis...");
        BufferedWriter buff = new BufferedWriter(new FileWriter(out));
        double start = System.currentTimeMillis();

        buff.write("term\tname\tsize\tis_true\tnoverlap\thg_pval\thg_fdr\tfej_pval\tfej_fdr\tks_stat\tks_pval\tks_fdr\tshortest_path_to_a_true");

        List<AnalysisEntry> entries = Collections.synchronizedList(new ArrayList<>());

        dag.getNodeMap().values().parallelStream().forEach(go -> {
            if (!(go.getGeneSymbols().size() >= min && go.getGeneSymbols().size() <= max))
            {
                return;
            }

            AnalysisEntry entry = new AnalysisEntry(go.getId());
            entry.setName(go.getAnnot());

            int n = smartIntersect(go.getGeneSymbols(), enrichedAll);
            entry.setSize(n);
            entry.setIs_true(go.isTrue());

            int k = smartIntersect(go.getGeneSymbols(), enrichedSigs);
            entry.setNoverlap(k);

            int N = smartIntersect(dag.getRoot().getGeneSymbols(), enrichedAll);
            int K = smartIntersect(dag.getRoot().getGeneSymbols(), enrichedSigs);

            double hgpval = calcHgPval(N, n, K, k);
            entry.setHg_pval(hgpval);

            double fejpval = calFejPval(N, n, K, k);
            entry.setFej_pval(fejpval);

            HashSet<String> bg = new HashSet<>(dag.getRoot().getGeneSymbols());
            HashSet<String> overlappingEnriched = intersect(go.getGeneSymbols(), enrichedAll);
            bg.removeAll(overlappingEnriched);

            double[] ksStats = calcKS(overlappingEnriched, bg);
            entry.setKs_stat(ksStats[0]);
            entry.setKs_pval(ksStats[1]);

            if (go.isTrue() || dag.getTrueGoEntries().isEmpty())
            {
                entry.setShortest_path_to_a_true("");
            }
            else
            {
                String path = go.getShortestPathToTrue(dag.getTrueGoIds(), dag);
                entry.setShortest_path_to_a_true(path);
            }

            entries.add(entry);
        });

        entries.sort(Comparator.comparing(AnalysisEntry::getTerm));

        List<Double> hgPvals = entries.stream().map(AnalysisEntry::getHg_pval).collect(Collectors.toList());
        List<Double> fejPvals = entries.stream().map(AnalysisEntry::getFej_pval).collect(Collectors.toList());
        List<Double> ksPvals = entries.stream().map(AnalysisEntry::getKs_pval).collect(Collectors.toList());

        List<Double> adjHgPvals = calculateBH_FDR(hgPvals);
        List<Double> adjFejPvals = calculateBH_FDR(fejPvals);
        List<Double> adjKsPvals = calculateBH_FDR(ksPvals);

        for (int i = 0; i < entries.size(); i++)
        {
            entries.get(i).setHg_fdr(adjHgPvals.get(i));
            entries.get(i).setFej_fdr(adjFejPvals.get(i));
            entries.get(i).setKs_fdr(adjKsPvals.get(i));

            buff.write("\n");
            buff.write(entries.get(i).toString());
        }

        buff.flush();
        logger.info(String.format("Time needed for Analysis: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
    }

    public HashSet<String> intersect(HashSet<String> A, HashSet<String> B)
    {
        HashSet<String> intersection = new HashSet<>(A);
        intersection.retainAll(B);
        return intersection;
    }

    public int smartIntersect(HashSet<String> A, HashSet<String> B)
    {
        int count = 0;
        for (String item : A)
        {
            if (B.contains(item)) count++;
        }
        return count;
    }


    public double calcHgPval(int N, int n, int K, int k)
    {
        HypergeometricDistribution hg = new HypergeometricDistribution(N, K, n);
        return hg.upperCumulativeProbability(k);
    }

    public double calFejPval(int N, int n, int K, int k)
    {
        HypergeometricDistribution hg = new HypergeometricDistribution(N - 1, K - 1, n - 1);
        return hg.upperCumulativeProbability(k - 1);
    }

    public ArrayList<Double> calculateBH_FDR(List<Double> pValues)
    {
        int m = pValues.size();
        ArrayList<Double> sortedPValues = new ArrayList<>(pValues);
        ArrayList<Double> bhAdjustedPValues = new ArrayList<>(Collections.nCopies(m, 1.0));

        HashMap<Double, ArrayList<Integer>> pValueToIndices = new HashMap<>();
        for (int i = 0; i < m; i++)
        {
            double pValue = pValues.get(i);
            pValueToIndices.computeIfAbsent(pValue, k -> new ArrayList<>()).add(i);
        }

        sortedPValues.sort(Double::compareTo);

        for (int i = 0; i < m; i++)
        {
            int rank = i + 1;
            double rawPValue = sortedPValues.get(i);
            double adjustedPValue = rawPValue * m / rank;
            bhAdjustedPValues.set(i, Math.min(adjustedPValue, 1.0));
        }

        for (int i = m - 2; i >= 0; i--)
        {
            if (bhAdjustedPValues.get(i) > bhAdjustedPValues.get(i + 1))
            {
                bhAdjustedPValues.set(i, bhAdjustedPValues.get(i + 1));
            }
        }

        ArrayList<Double> finalAdjustedPValues = new ArrayList<>(Collections.nCopies(m, 1.0));
        for (int i = 0; i < m; i++)
        {
            double pValue = sortedPValues.get(i);
            ArrayList<Integer> originalIndices = pValueToIndices.get(pValue);
            double adjustedPValue = bhAdjustedPValues.get(i);

            for (int originalIndex : originalIndices)
            {
                finalAdjustedPValues.set(originalIndex, adjustedPValue);
            }

            originalIndices.clear();
        }

        return finalAdjustedPValues;
    }

    public double[] calcKS(HashSet<String> in_set, HashSet<String> bg)
    {
        double[] inSetDist = new double[in_set.size()];
        double[] bgDist = new double[bg.size()];
        int i = 0;
        for (String geneId : in_set)
        {
            if (enrichedGeneMap.containsKey(geneId))
            {
                Gene curr = enrichedGeneMap.get(geneId);
                inSetDist[i] = curr.getFc();
                i++;
            }
        }

        i = 0;
        for (String geneId : bg)
        {
            if (enrichedGeneMap.containsKey(geneId))
            {
                Gene curr = enrichedGeneMap.get(geneId);
                bgDist[i] = curr.getFc();
                i++;
            }
        }
        double pval = this.ks.kolmogorovSmirnovTest(inSetDist, bgDist);
        double stat = this.ks.kolmogorovSmirnovStatistic(inSetDist, bgDist);
        return new double[]{stat, pval};
    }
}
