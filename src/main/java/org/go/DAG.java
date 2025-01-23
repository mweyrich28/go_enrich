package org.go;

import java.util.ArrayList;
import java.util.HashMap;

public class DAG
{

    private HashMap<String, GOEntry> nodeMap;
    private String type;


    // METHOD ADD NODE CHECK IF ALREADY IN H MAP
    public DAG(String type)
    {
        this.type = type;
        this.nodeMap = new HashMap<>();
    }

    public void insertNode(GOEntry go)
    {
        if (!nodeMap.containsKey(go.getId()))
        {
            nodeMap.put(go.getId(), go);
        }
    }

    @Override
    public String toString()
    {
        StringBuilder sb = new StringBuilder();
        for (String goId : nodeMap.keySet())
        {
            sb.append(goId);
            if (nodeMap.get(goId).getChildren().size() == 0)
            {
                sb.append("\t STOP");
            } else
            {
                sb.append("\tCHILDREN ->");
                sb.append("\t");
                ArrayList<String> ids = new ArrayList<>();
                for (GOEntry go : nodeMap.get(goId).getChildren())
                {
                    ids.add(go.getId());
                }
                sb.append(String.join(" | ", ids));
            }
            sb.append("\n\n");
        }

        return sb.toString();
    }

    public HashMap<String, GOEntry> getNodeMap()
    {
        return nodeMap;
    }
}
