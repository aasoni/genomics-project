package bf;
import java.util.*;
import java.io.*;
public class bloom_driver{

    public static class read
    {
		 /* read class:
		  * - name read1002
		 * - seq AGCT
		 * - kmer_len length of kmers used
		 * - num_kmers total
		 * - #kmers in +
		 * - #kmers in -
		 * - actual_pos T/F
		 * - called_pos T/F*/
    	
    	//constructor
 	   read( String namey, String seqy, int len )
 	   {
 		   name  	= namey;
 		   seq = seqy;
 		   kmer_len = len;
 	   }
 	   String name;
	   String seq;
	   int kmer_len;
	   int num_kmers = 0;
	   int num_pos_kmers = 0;
	   int num_neg_kmers = 0;
	   boolean actual_pos;
	   boolean called_pos;
	   /*
	   public void set_num_kmers(int num_kmersy) {
	       this.num_kmers = num_kmersy;
	   }
	   public int get_num_kmers() {
	       return this.num_kmers;
	   }
	   public void set_num_pos_kmers(int num_pos_kmersy) {
	       this.num_pos_kmers = num_pos_kmersy;
	   }
	   public int get_num_pos_kmers() {
	       return this.num_pos_kmers;
	   }
	   public void set_num_neg_kmers(int num_neg_kmersy) {
	       this.num_neg_kmers = num_neg_kmersy;
	   }
	   public int get_num_neg_kmers() {
	       return this.num_neg_kmers;
	   }
	   public void set_actual_pos(boolean actual_posy) {
	       this.actual_pos = actual_posy;
	   }
	   public boolean get_actual_pos() {
	       return this.actual_pos;
	   }
	   public void set_called_pos(boolean called_posy) {
	       this.called_pos = called_posy;
	   }
	   public boolean get_called_pos() {
	       return this.called_pos;
	   }
	   */
    }
	/**
	 * @param args
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws FileNotFoundException {
		/* 1) read in + and - enhancers and count 
		 * 2) save first 1/4 for testing,
		 * 3) use 3/4 of enhancers to break into kmers and train + and - bloom filters
		 * 4) run in test data 
		 * 5) report # correct and data
		 * */
		double false_positive_probability = 0.1;
		int kmer_size = 0;
		boolean normalize = false;
		int testcase = 0;
		if(args.length == 3)
		{
			kmer_size = Integer.parseInt(args[0]);
			testcase = Integer.parseInt(args[1]);
			normalize = Boolean.parseBoolean(args[2]);
		}
		else
		{
			System.out.println("usage: kmer_size(int) testcase(1-4) normalize(true/false)");
			System.exit(0);
		}
		Scanner infile1 = new Scanner( new FileReader( "enh_fb.fa" ) );
		Scanner infile2 = new Scanner( new FileReader( "nullseqsi_200_1.fa" ) );
		String title;
		String seq;
		int in_size = 5000;
		read[] pos_reads = new read[in_size];
		read[] neg_reads = new read[in_size];
		int num_pos_reads = 0;
		int num_neg_reads = 0;
		int total_pos_length = 0;
		int total_neg_length = 0;
		
		
		while(infile1.hasNext())
		{//read in positive naively
			title = infile1.nextLine();
			if(infile1.hasNext())
			{
				seq = infile1.nextLine();
				pos_reads[num_pos_reads++] = new read(title.substring(1), seq, kmer_size);
				total_pos_length += seq.length();
			}

			
		}
		while(infile2.hasNext())
		{//read in negative naively
			title = infile2.nextLine();
			if(infile2.hasNext())
			{
				seq = infile2.nextLine();
				neg_reads[num_neg_reads++] = new read(title.substring(1), seq, kmer_size);
				total_neg_length += seq.length();
			}
			
		}
		
		if(normalize)
		{ //ensure to use lower value
			if( num_neg_reads > num_pos_reads) num_neg_reads = num_pos_reads;
			else num_pos_reads = num_neg_reads; 
		}
		//test on 1/4, train on 3/4
		int[] ptest = new int[num_pos_reads / 4];
		int[] ntest = new int[num_neg_reads / 4];
		int[] ptrain = new int[num_pos_reads - num_pos_reads/4];
		int[] ntrain = new int[num_neg_reads - num_neg_reads/4];
		int pcount = 0; int ncount = 0; int pt = 0; int nt = 0;

		switch(testcase)
		{
			case 1: //test on first 1/4, train on rest
				for(int i = 0; i < num_pos_reads / 4; i++)
					ptest[pt++] = i;
				for(int i = 0; i < num_neg_reads / 4; i++)
					ntest[nt++] = i;
				for(int i = num_pos_reads/4; i < num_pos_reads; i++)
					ptrain[pcount++] = i;
				for(int i = num_neg_reads/4; i < num_neg_reads; i++)
					ntrain[ncount++] = i;
				break;
			case 2://test on second 1/4, train on rest
				for(int i = num_pos_reads / 4; i < (2 * num_pos_reads) / 4; i++)
					ptest[pt++] = i;
				for(int i = num_neg_reads / 4; i < (2 * num_neg_reads) / 4; i++)
					ntest[nt++] = i;
				for(int i = 0; i < num_pos_reads / 4; i++)
					ptrain[pcount++] = i;
				for(int j = (2 * num_pos_reads) / 4; j < num_pos_reads; j++)
					ptrain[pcount++] = j;
				for(int i = 0; i < num_neg_reads / 4; i++)
					ntrain[ncount++] = i;
				for(int j = (2 * num_neg_reads) / 4; j < num_neg_reads; j++)
					ntrain[ncount++] = j;
				break;
			case 3://test on third 1/4, train on rest
				for(int i = (2 * num_pos_reads) / 4; i < (3 * num_pos_reads) / 4; i++)
					ptest[pt++] = i;
				for(int i = (2 * num_neg_reads) / 4; i < (3 * num_neg_reads) / 4; i++)
					ntest[nt++] = i;
				for(int i = 0; i < (2 * num_pos_reads) / 4; i++)
					ptrain[pcount++] = i;
				for(int j = (3 * num_pos_reads) / 4; j < num_pos_reads; j++)
					ptrain[pcount++] = j;
				for(int i = 0; i < (2 * num_neg_reads) / 4; i++)
					ntrain[ncount++] = i;
				for(int j = (3 * num_neg_reads) / 4; j < num_neg_reads; j++)
					ntrain[ncount++] = j;
				break;
			case 4://test on fourth 1/4, train on rest
				for(int i = (3 * num_pos_reads) / 4; i < num_pos_reads; i++)
					ptest[pt++] = i;
				for(int i = (3 * num_neg_reads) / 4; i < num_neg_reads; i++)
					ntest[nt++] = i;
				for(int i = 0; i < (3 * num_pos_reads) / 4; i++)
					ptrain[pcount++] = i;
				for(int i = 0; i < (3 * num_neg_reads) / 4; i++)
					ntrain[ncount++] = i;
				break;
		}
		//real training length is only 3/4 so calculate
		int use_pos_length = 0; int use_neg_length = 0;
		for(int i = 0; i < ptrain.length; i++) 
		{
			use_pos_length += pos_reads[ptrain[i]].seq.length();
		}

		
		for(int i = 0; i < ntrain.length; i++) 
		{
			use_neg_length += neg_reads[ntrain[i]].seq.length();
		}
		
		BloomFilter<String> bloom_filter_pos = new BloomFilter<String>(false_positive_probability, use_pos_length);
		BloomFilter<String> bloom_filter_neg = new BloomFilter<String>(false_positive_probability, use_neg_length);
		
		String kmer = "";
		//training pos and neg filters
		for(int i = 0; i < ptrain.length; i++) 
		{//for training portion of the reads
			for(int j = 0; j <= pos_reads[ptrain[i]].seq.length() - kmer_size; j++) 
			{//for each character make kmer of kmer_size
				kmer = pos_reads[ptrain[i]].seq.substring(j, j + kmer_size);
				bloom_filter_pos.add(kmer);	
				
			}
			pos_reads[ptrain[i]].num_kmers = pos_reads[ptrain[i]].seq.length() - kmer_size + 1;
			pos_reads[ptrain[i]].actual_pos = true;
		}
		for(int i = 0; i < ntrain.length; i++) 
		{//for training portion of the reads
			for(int j = 0; j <= neg_reads[ntrain[i]].seq.length() - kmer_size; j++) 
			{//for each character make kmer of kmer_size
				kmer = neg_reads[ntrain[i]].seq.substring(j, j + kmer_size);
				bloom_filter_neg.add(kmer);	
			}
			neg_reads[ntrain[i]].num_kmers = neg_reads[ntrain[i]].seq.length() - kmer_size + 1;
			neg_reads[ntrain[i]].actual_pos = false;
		}

		//test on 1/4 of saved data - the test set
		int num_pos_right = 0; //number of positive reads called positive
		int num_neg_right = 0; //number of negative reads called negative
		double pos_neg_diff = 0; //absolute difference between # kmers called pos and neg
		
		for(int i = 0; i < ptest.length; i++) 
		{//for test portion of the reads
			for(int j = 0; j <= pos_reads[ptest[i]].seq.length()-kmer_size; j++) 
			{//for each character make kmer of kmer_size
				kmer = pos_reads[ptest[i]].seq.substring(j, j+kmer_size);
				if(bloom_filter_pos.contains(kmer))
				{
					pos_reads[ptest[i]].num_pos_kmers++;
				}
				if(bloom_filter_neg.contains(kmer))
				{
					pos_reads[ptest[i]].num_neg_kmers++;	
				}
			}
			pos_reads[ptest[i]].num_kmers = pos_reads[ptest[i]].seq.length()-kmer_size;
			pos_reads[ptest[i]].actual_pos = true;
			if( pos_reads[ptest[i]].num_pos_kmers - pos_reads[ptest[i]].num_neg_kmers > 0)
			{
				pos_neg_diff += pos_reads[ptest[i]].num_pos_kmers - pos_reads[ptest[i]].num_neg_kmers;
				pos_reads[ptest[i]].called_pos = true;
				num_pos_right += 1;
			}
			else //more negative so call negative
			{
				pos_neg_diff += pos_reads[ptest[i]].num_neg_kmers - pos_reads[ptest[i]].num_pos_kmers;
				pos_reads[ptest[i]].called_pos = false;
			}
		}
		//negative tests
		for(int i = 0 ; i < ntest.length; i++) 
		{//for first 1/4 of the reads
			for(int j = 0; j <= neg_reads[ntest[i]].seq.length()-kmer_size; j++) 
			{//for each character make kmer of kmer_size
				kmer = neg_reads[ntest[i]].seq.substring(j, j+kmer_size);
				if(bloom_filter_pos.contains(kmer))
				{
					neg_reads[ntest[i]].num_pos_kmers++;
				}
				if(bloom_filter_neg.contains(kmer))
				{
					neg_reads[ntest[i]].num_neg_kmers++;
				}
			}
			neg_reads[ntest[i]].num_kmers = neg_reads[ntest[i]].seq.length()-kmer_size;
			neg_reads[ntest[i]].actual_pos = false;
			if( neg_reads[ntest[i]].num_pos_kmers - neg_reads[ntest[i]].num_neg_kmers > 0)
			{
				pos_neg_diff += neg_reads[ntest[i]].num_pos_kmers - neg_reads[ntest[i]].num_neg_kmers;
				neg_reads[ntest[i]].called_pos = true;
			}
			else //more negative so call negative
			{
				pos_neg_diff += neg_reads[ntest[i]].num_neg_kmers - neg_reads[ntest[i]].num_pos_kmers;
				neg_reads[ntest[i]].called_pos = false;
				num_neg_right += 1;
			}
		}
		double average_pos_neg_difference = pos_neg_diff / (ptest.length + ntest.length);
		System.out.println("kmer_size: " + kmer_size);
		System.out.println("#positive right: " + num_pos_right);
		System.out.println("#positive tests: " + ptest.length);
		System.out.println("#negative right: " + num_neg_right);
		System.out.println("#negtaive tests: " + ntest.length);
		System.out.println("average_pos_neg_difference: " + average_pos_neg_difference);
		double average_kmers = (total_pos_length + total_neg_length) / (num_pos_reads + num_neg_reads);
		System.out.println("average_kmers per read: " + average_kmers);

	}

}
