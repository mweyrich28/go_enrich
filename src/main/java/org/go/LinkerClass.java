package org.go;

public class LinkerClass
{
    private GOEntry go;
    private int k;

    public LinkerClass(GOEntry go, int k)
    {
        this.go = go;
        this.k = k;
    }

    public GOEntry getGo()
    {
        return go;
    }

    public void setGo(GOEntry go)
    {
        this.go = go;
    }

    public int getK()
    {
        return k;
    }

    public void setK(int k)
    {
        this.k = k;
    }
}
