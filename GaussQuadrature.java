import java.util.Arrays;

public class GaussQuadrature {

    public static void Gauss()
    {
        System.out.println("Guass Quadrature 5-Point");
        /***
         * setting up values for nodes1 and weights
         */
        double[] nodes1 = new double[5];
        double[] nodes2 = new double[5];
        double[] nodes3 = new double[5];
        double[] weights = new double[5];
        //assigning values to nodes1
        double a = Math.sqrt(  (1.0/9.0)* (5- (2* (Math.sqrt(10.0/7.0)) )) );
        double b = Math.sqrt(  (1.0/9.0)* (5 + (2* (Math.sqrt(10.0/7.0)) )) );
        nodes1[0]= -a;
        nodes1[1] = -b;
        nodes1[2] = 0;
        nodes1[3] = a;
        nodes1[4] = b;
        //assigning same values to nodes2
        for(int i = 0; i <nodes2.length; i++)
        {
            nodes2[i] = nodes1[i];
        }
        //assigning same values to nodes3
        for(int i = 0; i <nodes3.length; i++)
        {
            nodes3[i] = nodes1[i];
        }
        //assigning values to weights
        double c = 0.3*( (-0.7 + (5.0 * Math.sqrt(0.7))) / (-2.0 + (5.0 * Math.sqrt(0.7) ) ) );
        double d = 0.3*( (0.7 + (5.0 * Math.sqrt(0.7))) / (2.0 + (5.0 * Math.sqrt(0.7) ) ) );
        weights[0] = c;
        weights[1] = d;
        weights[2] = (128.0/225.0);
        weights[3]=c;
        weights[4]=d;
        System.out.println("Weights: " + Arrays.toString(weights));
        System.out.println("Nodes: " + Arrays.toString(nodes1));
        /***
         * We now have Ai and xi arrays, with xi arrays (node arrays) for 3 functions. Now lets pass xi through
         * each function (nodes1 through f1, nodes2 through f2, nodes3 through f3)
         */
        Functions.f1(nodes1);
        Functions.f2(nodes2);
        Functions.f3(nodes3);
        /***
         * Now lets hold the values of nodes(1,2,3) * weights. Remember nodes(1,2,3) now hold the nodes passed
         * through functions, so it now holds f(xi) for all three functions in respective nodes(1,2,3) array.
         */
        for(int i = 0; i<nodes1.length; i++)
        {
            nodes1[i] = nodes1[i] * weights[i];
        }
        for(int i = 0; i<nodes2.length; i++)
        {
            nodes2[i] = nodes2[i] * weights[i];
        }
        for(int i = 0; i<nodes3.length; i++)
        {
            nodes3[i] = nodes3[i] * weights[i];
        }
        /***
         * All thats left is to sum up the Ai * f(xi) terms and we get our approximation
         */
        double func1Appx = 0;
        double func2Appx = 0;
        double func3Appx = 0;
        for(int i = 0; i<nodes1.length; i++)
        {
            func1Appx = func1Appx + nodes1[i];
        }
        for(int i = 0; i<nodes2.length; i++)
        {
            func2Appx = func2Appx + nodes2[i];
        }
        for(int i = 0; i<nodes3.length; i++)
        {
            func3Appx = func3Appx + nodes3[i];
        }
        //True values of functions, TV = True Value
        double f1TV = 0.58385316345285761300;
        double f2TV = 3.8856180831641267317;
        double f3TV = 1.3250027473578644309;

        System.out.println("Appx of f1: " + func1Appx + " True Value: " +
                f1TV + " Error: " + (f1TV - func1Appx));
        System.out.println("Appx of f2: " + func2Appx + " True Value: " +
                f2TV + " Error: " + (f2TV - func2Appx));
        System.out.println("Appx of f3: " + func3Appx + " True Value: " +
                f3TV + " Error: " + (f3TV - func3Appx));
        System.out.println("--------------------------------------------------------" +
                "-------------------------------------------------------------");
    }
}
