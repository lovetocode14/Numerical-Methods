public class RangeKutta {
    //run keeps track of run # for header
    static int run =1;
    /***
     * all coefficients used
     */
    static double alpha1 = .4;
    static double alpha2 = .8;
    static double alpha3 = 1.0;
    static double c1 = (-1 + (2*alpha2))/(12*(-1+alpha1)*alpha1*(alpha1-alpha2));
    static double c2 = (1-(2*alpha1))/(12*(alpha1-alpha2)*(-1+alpha2)*alpha2);
    static double c3 = (3-(4*alpha2) + (alpha1*(-4+(6*alpha2))))/(12*(-1+alpha1)*(-1+alpha2));
    static double c0 = 1-c1-c2-c3;
    static double beta10 = alpha1;
    static double beta21 = ((alpha1-alpha2)*alpha2)/(alpha1*(-2+(4*alpha1)));
    static double beta20 = alpha2 - beta21;
    static double beta31 = ((-1 + alpha1)*(-2+alpha1 +((5 - (4*alpha2))*alpha2)))/(2*alpha1*(alpha1-alpha2)*(3-(4*alpha2)+(alpha1*(-4+(6*alpha2)))));
    static double beta32 = ((-1+alpha1)*(-1+(2*alpha1))*(-1+alpha2))/((alpha1-alpha2)*alpha2*(3-(4*alpha2)+(alpha1*(-4+(6*alpha2)))));
    static double beta30 = alpha3 - beta31 - beta32;
    /***
     * function 3x^2 * cos(5t)
     * @param t
     * @param x
     * @return
     */
    private static double f(double t, double x)
    {
        return (3.0*Math.pow(x, 2))* Math.cos(5*t);
    }

    /* Though about defining a setCoefficients method here,
        to encapsulate the coefficients, and so that an end
        user could set their own, but using the dot operator
        or getters seemed tedious.
     */
    private static void rangeKutta(double h)
    {
        /***
         * steps = number of steps (range/h)
         * x0, t0, intial values
         * x1, x2, x3, x4, four step x values
         * trueAnswer for error calculcation
         * error variable
         */
        double steps = 2.0/h;
            double x0 = 1.0;
            double t0 = 0.0;
            double x1, x2, x3, x4;
            double trueAnswer = 0.753913186469598763502963347;
            double error;


            /*Simple table setup with inital values before they are
                changed, a bit rudimentary.
             */
            System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
            System.out.println("Run " + run);
            System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
            System.out.println("t + h               " + "xn           " + "(h=" + h + ")");
            System.out.println("-----------------------------------------");
            System.out.println(t0 + "                 " + x0);

            //Main logic, implementing the four stage RK4 summation formulas
            for(int j = 0; j< steps; j++)
            {
                x1 = x0 + h*(beta10*f(t0, x0));
                x2 = x0 + h*((beta20*f(t0, x0)) + (beta21*f(t0 + alpha1*h, x1)));
                x3 = x0 + h*((beta30*f(t0, x0)) + (beta31*f(t0 + alpha1*h, x1)) + (beta32*f(t0 + alpha2*h, x2)));
                x4 = x0 + h*((c0*f(t0, x0)) + (c1*f(t0 + alpha1*h, x1)) + (c2*f(t0 + alpha2*h, x2)) + (c3*f(t0 + alpha3*h, x3)));

                //Set t and x for next step
                t0 = t0+h;
                x0 = x4;

                //Display t and x values in table
                System.out.format("%1.7f", t0);
                System.out.println("           " + x0);
            }
        //calculate error and display
        error = Math.abs(trueAnswer - x0);
        System.out.println("-----------------------------------------");
        System.out.println("Error               " +error);
        //run counter for header
        run++;
    }

    public static void main(String[] args)
    {
        rangeKutta(0.1);
        rangeKutta(.05);
        rangeKutta(0.025);
        rangeKutta(0.0125);
        rangeKutta(0.00625);
        rangeKutta(0.003125);
        rangeKutta(0.0015625 );
    }
}
