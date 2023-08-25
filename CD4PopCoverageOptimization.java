import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Scanner;

public class PopCoverageOptimization {

    public static void main(String[] args)
    {
//         String data = "HLA-B*40:01	AAEIRASANL\r\n" //THIS IS WHERE YOU PASTE YOUR DATA FROM EXCEL. EXAMPLE PROVIDED
//         		+ "HLA-A*02:03	VLNDIFSRL\r\n";
        String data = args[0];

        String[] rowArray = data.split("\\r?\\n");
        String[] alleles = new String[rowArray.length];
        String[] epitopes = new String[rowArray.length];

        for(int i = 0; i < rowArray.length; i++)
        {
            alleles[i] = rowArray[i].substring(0,14);
            epitopes[i] = rowArray[i].substring(rowArray[i].lastIndexOf("\t")+1);
        }

        String[] epitopesNR = new HashSet<String>(Arrays.asList(epitopes)).toArray(new String[0]);
        String[][] alleleList = new String[epitopesNR.length][0];

        for(int i = 0; i < epitopesNR.length; i++)
        {
            String tempEpi = epitopesNR[i];

            for(int w = 0; w < epitopes.length; w++)
            {
                if(epitopes[w].equals(tempEpi))
                {
                    String[] tempArray = new String[alleleList[i].length+1];
                    for(int k = 0; k < alleleList[i].length; k++)
                    {
                        tempArray[k] = alleleList[i][k];
                    }
                    tempArray[tempArray.length-1] = alleles[w];
                    alleleList[i] = tempArray;
                }
            }
        }

        ArrayList<String> optimEpi = new ArrayList<String>();
        HashSet<String> allelesNRARRAY = new HashSet<String>(Arrays.asList(alleles));
        ArrayList<String> allelesNR = new ArrayList<String>();
        for(String s:allelesNRARRAY)
        {
            if(!allelesNR.contains(s))
            {
                allelesNR.add(s);
            }
        }
        ArrayList <String[]> optimAllele = new ArrayList<String[]>();

        String searchAllele;
        ArrayList<Integer> matchIndex = new ArrayList<Integer>();
        int index = 0;
        int temp = 0;

        while(allelesNR.size() > 0)
        {
            searchAllele = allelesNR.get(0);
            for(int i = 0; i < alleleList.length; i++)
            {
                for(int w = 0; w < alleleList[i].length; w++)
                {
                    if(alleleList[i][w].equals(searchAllele))
                    {
                        matchIndex.add(i);
                    }
                }
            }
            for(int i = 0; i < matchIndex.size(); i++)
            {
                if(alleleList[matchIndex.get(i)].length > temp)
                {
                    temp = alleleList[matchIndex.get(i)].length;
                    index = matchIndex.get(i);
                }
            }
            optimEpi.add(epitopesNR[index]);
            optimAllele.add(alleleList[index]);

			/*for(String[] s1:optimAllele)
			{
				for(String s2:s1)
				{
					System.out.print(s2 + " ");
				}
				System.out.println();
			}*/

            for(String s:alleleList[index])
            {
                allelesNR.remove(s);
            }
            matchIndex.clear();
            index = 0;
            temp = 0;
        }

        for(int i = 0; i < optimEpi.size(); i++)
        {
            System.out.print(optimEpi.get(i) + " ");
            for(int w = 0; w < optimAllele.get(i).length-1; w++)
            {
                System.out.print(optimAllele.get(i)[w] + ",");
            }
            System.out.print(optimAllele.get(i)[optimAllele.get(i).length-1]);
            System.out.println();
        }
    }

}