package org.go.util;

import org.go.DAG;
import org.go.GOEntry;

import java.io.*;
import java.util.ArrayList;
import java.util.zip.GZIPInputStream;

public class FileUtils
{

    public static DAG parseObo(String pathToObo, String type) throws IOException
    {
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
                for (String parentString : parents) // update its parents
                {
                    String[] comp = parentString.split(" ! ");
                    String parentId = comp[0].substring(6);
                    if (dag.getNodeMap().containsKey(parentId)) // check if parent already exists
                    {
                        dag.getNodeMap().get(parentId).addChild(currentGO);
                    }
                    else // create parent and then add currentGO
                    {
                        String parentAnnot = comp[1];
                        GOEntry newParent = new GOEntry(parentId, parentAnnot);
                        dag.insertNode(newParent);
                        newParent.addChild(currentGO);
                    }
                }
                // clear up vars
                parents.clear();
                currNodeId = null;
                namespace = null;
                desc = null;
            }
            else if (line.charAt(0) == 'i' && line.charAt(3) == 'o') // ignore entire entry
            {
                currNodeId = null;
                namespace = null;
                parents.clear();
            }
            else if (line.charAt(0) == 'i' && line.charAt(1) == 'd') // hit new entry id
            {
                currNodeId = line.split(" ")[1];
            }
            else if (line.charAt(0) == 'i' && line.charAt(3) == 'a') // hit new parent of entry id
            {
                parents.add(line);
            }
            else if (line.charAt(0) == 'n' && line.charAt(5) == 'p') // hit new parent of entry id
            {
                namespace = line.split(" ")[1];
            }
            else if (line.charAt(0) == 'n' && line.charAt(4) == ':') // hit desc of node
            {
                desc = line.split(": ")[1];
            }
        }
        System.out.println(dag);
        return dag;
    }

    public static void parseGAF(String pathToGaf) throws IOException
    {
        FileInputStream fileInputStream = new FileInputStream(pathToGaf);
        GZIPInputStream gzipInputStream = new GZIPInputStream(fileInputStream);
        InputStreamReader inputStreamReader = new InputStreamReader(gzipInputStream);
        BufferedReader bufferedReader = new BufferedReader(inputStreamReader);
        String line;
        while ((line = bufferedReader.readLine()) != null)
        {
            System.out.println(line);
        }
    }

    public static void parseENSG(String pathToENSG) throws IOException
    {

    }

}
