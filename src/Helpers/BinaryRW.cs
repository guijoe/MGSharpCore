using System;
using System.IO;
using System.ComponentModel;
using System.Runtime.Serialization.Formatters.Binary;

[Serializable()]
public class Loan : INotifyPropertyChanged
{
    public double LoanAmount { get; set; }
    public double InterestRate { get; set; }

    [field: NonSerialized()]
    public DateTime TimeLastLoaded { get; set; }

    const string FileName = "SavedLoan.bin";

    public int Term { get; set; }

    private string customer;
    public string Customer
    {
        get { return customer; }
        set
        {
            customer = value;
            //PropertyChanged?.Invoke(this,
            // new PropertyChangedEventArgs(nameof(Customer)));
        }
    }

    [field: NonSerialized()]
    public event System.ComponentModel.PropertyChangedEventHandler PropertyChanged;

    public Loan(double loanAmount,
                double interestRate,
                int term,
                string customer)
    {
        this.LoanAmount = loanAmount;
        this.InterestRate = interestRate;
        this.Term = term;
        this.customer = customer;
    }

    public static void Main(string[] args)
    {
        Loan TestLoan = new Loan(10000.0, 0.075, 36, "Neil Black");

        if (File.Exists(FileName))
        {
            Console.WriteLine("Reading saved file");
            Stream openFileStream = File.OpenRead(FileName);
            BinaryFormatter deserializer = new BinaryFormatter();
            TestLoan = (Loan)deserializer.Deserialize(openFileStream);
            TestLoan.TimeLastLoaded = DateTime.Now;
            openFileStream.Close();
        }

        //TestLoan.PropertyChanged += (_, __) => Console.WriteLine($"New customer value: {TestLoan.Customer}");

        TestLoan.Customer = "Henry Clay";
        Console.WriteLine(TestLoan.InterestRate);
        TestLoan.InterestRate = 7.1;
        Console.WriteLine(TestLoan.InterestRate);

        Stream SaveFileStream = File.Create(FileName);
        BinaryFormatter serializer = new BinaryFormatter();
        serializer.Serialize(SaveFileStream, TestLoan);
        SaveFileStream.Close();

        Console.ReadLine();
    }
}