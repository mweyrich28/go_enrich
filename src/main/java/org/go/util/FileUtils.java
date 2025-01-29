package org.go.util;

import org.go.DAG;
import org.go.GOEntry;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

import org.go.Gene;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class FileUtils
{

    private static final Logger logger = LoggerFactory.getLogger(FileUtils.class);

    public static DAG parseOBO(String pathToObo, String type) throws IOException
    {

        logger.info("Starting to parse obo...");
        double start = System.currentTimeMillis();
        BufferedReader bufferedReader = new BufferedReader(new FileReader(new File(pathToObo)));
        String line;
        String currNodeId = null;
        String namespace = null;
        String desc = null;
        ArrayList<String> parents = new ArrayList<>();
        DAG dag = new DAG(type);
        while ((line = bufferedReader.readLine()) != null)
        {
            if (line.isEmpty())
            {
                continue;
            }
            // we hit new entry
            if (line.charAt(0) == '[' && currNodeId != null && namespace.equals(type))
            {
                GOEntry currentGO;
                if (dag.getNodeMap().containsKey(currNodeId)) // check if node is already in dag
                {
                    currentGO = dag.getNodeMap().get(currNodeId);
                }
                else
                {
                    currentGO = new GOEntry(currNodeId, desc);
                    dag.insertNode(currentGO);
                }

                // check if we got root
                if (desc.equals(type))
                {
                    dag.setRoot(currentGO);
                }

                for (String parentString : parents) // update its parents
                {
                    String[] comp = parentString.split(" ! ");
                    String parentId = comp[0].substring(6);
                    if (dag.getNodeMap().containsKey(parentId)) // check if parent already exists
                    {
                        dag.getNodeMap().get(parentId).addChild(currentGO);

                        // debug: add parents
                        currentGO.addParent(dag.getNodeMap().get(parentId));
                    }
                    else // create parent and then add currentGO
                    {
                        String parentAnnot = comp[1];
                        GOEntry newParent = new GOEntry(parentId, parentAnnot);
                        dag.insertNode(newParent);
                        newParent.addChild(currentGO);

                        // debug: add parents
                        currentGO.addParent(newParent);
                    }
                }
                // clear up vars
                parents.clear();
                currNodeId = null;
                namespace = null;
                desc = null;
            }
            // reset vars if above condition is not met
            else if (line.startsWith("["))
            {
                parents.clear();
                currNodeId = null;
                namespace = null;
                desc = null;
            }
            // ignore entire entry
            else if (line.charAt(0) == 'i' && line.charAt(3) == 'o')
            {
                currNodeId = null;
                namespace = null;
                desc = null;
                parents.clear();
            }
            // hit new entry id
            else if (line.charAt(0) == 'i' && line.charAt(1) == 'd')
            {
                currNodeId = line.split(" ")[1];
            }
            // hit new parent of entry id
            else if (line.charAt(0) == 'i' && line.charAt(3) == 'a')
            {
                parents.add(line);
            }
            // hit new parent of entry id
            else if (line.charAt(0) == 'n' && line.charAt(5) == 'p')
            {
                namespace = line.split(" ")[1];
            }
            // hit desc of node
            else if (line.charAt(0) == 'n' && line.charAt(4) == ':')
            {
                desc = line.split(": ")[1];
            }
        }

        logger.info(String.format("Time needed for parsing obo: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
        logger.info("Created DAG contains " + dag.getNodeMap().size() + " GO entries");
        return dag;
    }

    public static void parseGAF(String pathToGaf, DAG dag) throws IOException
    {
        logger.info("Starting to parse gaf...");
        double start = System.currentTimeMillis();
        FileInputStream fileInputStream = new FileInputStream(pathToGaf);
        GZIPInputStream gzipInputStream = new GZIPInputStream(fileInputStream);
        InputStreamReader inputStreamReader = new InputStreamReader(gzipInputStream);
        BufferedReader bufferedReader = new BufferedReader(inputStreamReader);
        String line;
        while ((line = bufferedReader.readLine()) != null)
        {
            if (line.charAt(0) == '!') // skip comments
            {
                continue;
            }
            String[] comp = line.split("\t");
            if(!comp[3].isEmpty()) {
                continue;
            }
            String geneSymbol = comp[2];
            String goName = comp[4];
            if (dag.getNodeMap().containsKey(goName))
            {
                dag.getNodeMap().get(goName).addSymbol(geneSymbol);
            }
        }
        logger.info(String.format("Time needed for parsing gaf: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
    }

    public static void parseENSEMBL(String pathToEnsembl, DAG dag) throws IOException
    {

        logger.info("Starting to parse ensemble...");
        double start = System.currentTimeMillis();
        BufferedReader bufferedReader = new BufferedReader(new FileReader(pathToEnsembl));
        String line = null;
        boolean skippedFirstLine = false;
        while ((line = bufferedReader.readLine()) != null)
        {
            if (!skippedFirstLine)
            {
                skippedFirstLine = true;
                continue;
            }

            StringBuilder sb = new StringBuilder();
            boolean seenFirstTab = false;
            boolean seenSecondTab = false;
            String currGeneSymbol = null;
            for (int i = 0; i < line.length(); i++)
            {
                if (line.charAt(i) == '\t' && !seenFirstTab)
                {
                    seenFirstTab = true;
                }
                else if (seenFirstTab && line.charAt(i) == '\t')
                {
                    seenSecondTab = true;
                    // check if we even gave a gene symbol
                    if (sb.isEmpty())
                    {
                        continue;
                    }
                    currGeneSymbol = sb.toString();
                    sb.setLength(0);
                }
                else if (seenFirstTab && !seenSecondTab)
                {
                    sb.append(line.charAt(i));
                }
                else if (seenSecondTab && line.charAt(i) != '|')
                {
                    sb.append(line.charAt(i));
                }
                if (seenSecondTab && (line.charAt(i) == '|' || i == line.length() - 1))
                {
                    String currentGo = sb.toString();
                    if (dag.getNodeMap().containsKey(currentGo) && currGeneSymbol != null)
                    {
                        dag.getNodeMap().get(currentGo).addSymbol(currGeneSymbol);
                    }
                    sb.setLength(0);
                }
            }
        }
        logger.info(String.format("Time needed for parsing ensmble: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
    }

    public static HashMap<String, Gene> parseEnrichedGenes(String pathToEnrichedGenes, DAG dag) throws IOException
    {
        logger.info("Starting to parse enriched genes...");
        double start = System.currentTimeMillis();
        HashMap<String, Gene> enrichedGenes = new HashMap<>();
        BufferedReader bufferedReader = new BufferedReader(new FileReader(pathToEnrichedGenes));

        String line;
        boolean skippedHeader = false;
        while ((line = bufferedReader.readLine()) != null)
        {
            if (line.charAt(0) == '#') // handle ground truth by updating go nodes
            {
                String trueGoId = line.substring(1);
                if (dag.getNodeMap().containsKey(trueGoId))
                {
                    dag.getNodeMap().get(trueGoId).setTrue();
                }
                else
                {
                    logger.debug("True GO Entry " + trueGoId + " not in DAG!");
                }
                continue;
            }
            else if (!skippedHeader) // skip header
            {
                skippedHeader = true;
                continue;
            }

            // parse remaining genes
            String[] comp = line.split("\t");
            if (comp[2].charAt(0) == 'f')
            {
                enrichedGenes.put(comp[0], new Gene(comp[0], Double.parseDouble(comp[1]), false));
            }
            else
            {
                enrichedGenes.put(comp[0], new Gene(comp[0], Double.parseDouble(comp[1]), true));
            }

        }

        logger.info(String.format("Time needed for parsing enriched genes: %s seconds", (System.currentTimeMillis() - start) / 1000.0));
        return enrichedGenes;
    }
}
