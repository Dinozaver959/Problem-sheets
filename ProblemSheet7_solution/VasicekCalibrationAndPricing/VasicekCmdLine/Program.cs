using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace VasicekCmdLine
{
    class Program
    {
        static void Main(string[] args)
        {
            Examples.TestVasicekPricing();
            Examples.TestVasicekCalibration();
            Console.ReadKey();
        }
    }
}
