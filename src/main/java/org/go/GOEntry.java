package org.go;

import java.util.ArrayList;
import java.util.HashSet;

public class GOEntry
{
    private String id;
    private String annot;
    private boolean isTrue;
    private ArrayList<GOEntry> parents;
    private ArrayList<GOEntry> children;
    private HashSet<String> geneSymbols;
    private int depth;

    public GOEntry(String id, String annot)
    {
        this.id = id;
        this.annot = annot;
        this.geneSymbols = new HashSet<>();
        this.parents = new ArrayList<>();
        this.children = new ArrayList<>();
    }

    public void setAnnot(String annot)
    {
        if (this.annot == null)
        {
            this.annot = annot;
        }
    }

    public String getId()
    {
        return id;
    }

    public void addSymbol(String geneId)
    {
        this.geneSymbols.add(geneId);
    }

    public void addParent(GOEntry parent)
    {
        this.parents.add(parent);
    }

    public void addChild(GOEntry child)
    {
        this.children.add(child);
    }

    public ArrayList<GOEntry> getChildren()
    {
        return children;
    }

    @Override
    public String toString()
    {
        return this.depth + "\t" + this.id + "\t" + this.annot + "\t" + children.size();
    }

    public void setTrue()
    {
        isTrue = true;
    }

    public void inheritGeneSymbols()
    {
        for (GOEntry child: children)
        {
            child.inheritGeneSymbols();
            this.getGeneSymbols().addAll(child.geneSymbols);
        }
    }

    public HashSet<String> getGeneSymbols()
    {
        return geneSymbols;
    }

    public void setDepth(int depth)
    {
        this.depth = depth;
    }
    public void passDepth()
    {
        int newDepth = this.depth;
        newDepth++;
        for (GOEntry child : children)
        {
            child.setDepth(newDepth);
            child.passDepth();
        }
    }

    public String getAnnot()
    {
        return annot;
    }

    public boolean isTrue()
    {
        return isTrue;
    }
}
