#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <math.h>


using namespace std;

static vector<int> private_key;
static vector<vector<int> > U;
static vector<vector<int> > HybridMatrix;
static int nblock_bio;
static int nblock_pri;
static vector<vector<int> > Matrix;

static vector<int> HYBRID_CODE;

//generation of hybrid_code and matrix and all
void hybrid_generation()
{
	for (int i = 0; i < HybridMatrix.size(); i++)
	{
		for (int j = 0; j < HybridMatrix[i].size(); j++)
		{
			HYBRID_CODE.push_back(HybridMatrix[i][j]);
		}
	}

	fstream hyb("hybrid_code.txt", ios::out);
	for (int i = 0; i < HYBRID_CODE.size(); i++)
	{
		hyb << HYBRID_CODE[i] << endl;
	}
	hyb.close();

}

void multiplication(vector<vector<int> > permMatrix)
{
	HybridMatrix.resize(U.size()); //declaring no. of rows to be same
	int i, j, k;
	for (i = 0; i < U.size(); i++)
	{
		for (j = 0; j < U.size(); j++)
		{
			HybridMatrix[i].push_back(0);
			for (k = 0; k < U.size(); k++)
				HybridMatrix[i][j] += U[i][k] * permMatrix[k][j];
		}
	}
}


void permutation_and_multiplication()
{
    ofstream perm("perm.txt");
	//6 permutation each 3x3 size
	int p0[3][3] = { {1,0,0},{0,1,0},{0,0,1} };
	int p1[3][3] = { {1,0,0},{0,0,1},{0,1,0} };
	int p2[3][3] = { {0,1,0},{1,0,0},{0,0,1} };
	int p3[3][3] = { {0,0,1},{1,0,0},{0,1,0} };
	int p4[3][3] = { {0,1,0},{0,0,1},{1,0,0} };
	int p5[3][3] = { {0,0,1},{0,1,0},{1,0,0} };

	int row_covered = 0;
	vector<vector<int> > permMatrix(U.size(), vector<int>(U.size(), 0));  // declare a matrix of same dimension as U, and initialize all elements as 0

	for (int i = 0; i < private_key.size(); i++)
	{
		int perm_select = private_key[i] % 6; //to select one of the above three permutations
		if (row_covered < U.size())
		{
			int start = 3 * i;
			int endd = (3 * i) + 3;
			for (int j = start; j < endd; j++)
			{
				for (int k = start; k < endd; k++)
				{
					if (perm_select == 0)
						permMatrix[j][k] = p0[j - (3 * i)][k - (3 * i)];
					else if (perm_select == 1)
						permMatrix[j][k] = p1[j - (3 * i)][k - (3 * i)];
					else if (perm_select == 2)
						permMatrix[j][k] = p2[j - (3 * i)][k - (3 * i)];
					else if (perm_select == 3)
						permMatrix[j][k] = p3[j - (3 * i)][k - (3 * i)];
					else if (perm_select == 4)
						permMatrix[j][k] = p4[j - (3 * i)][k - (3 * i)];
					else if (perm_select == 5)
						permMatrix[j][k] = p5[j - (3 * i)][k - (3 * i)];

				}
			}
			perm << perm_select << endl; // writing in file, which permutation is selected
			row_covered += 3; // 3 rows covered
		}
		else
			break;
	}

	//writing permutation matrix int a file
	perm.close();

	multiplication(permMatrix);
}


void addExtra(int &last_row_filled, int col)
{
    int *point=private_key.data();
    int last_row_col=Matrix[last_row_filled].size();
    int count=0;

    if(last_row_col==col) //if all columns are filled
    {
        last_row_filled++;
    }
    if(last_row_col!=col)
    {
        for(int i=last_row_col;i<col;i++)
        {
            if(count<private_key.size())
            {
                //int index=(*point++)%(col);
                int index=*point++; //modules is taken so as to not make vector cross boundary
                count++;
                if(i%2==0)
                    Matrix[last_row_filled].push_back(Matrix[0][index]);
                else
                    Matrix[last_row_filled].push_back(Matrix[index][0]);
            }
            else
            {
                count=0;
                point=private_key.data();
                continue;
            }
        }
        last_row_filled++;
    }
}


vector<vector<int> > matrix_generation(int first_ele, int second_ele, vector<int> first, vector<int> second, int block1, int block2)
{
	int matrix_dimen = ceil(sqrt(first.size() + second.size())); //dimension of matrix
	int tmp_blk_one = block1; // to keep in check the blocks used of first
	int tmp_blk_two = block2; // to keep in check the blocks used of second

	Matrix.resize(matrix_dimen);
	int* pointFirst = first.data(); //points to the first element of first vector
	int* pointSecond = second.data(); //points to the second element of first vector

	int cell_fill = 1; //cells to be filled
	int row_to_fill = 0; //rows to be filled
	int count = 0; //keep check on private_key

	for (int i = 0; i < private_key.size(); i++) //run full private key vector
	{
		if (count < private_key.size())
		{
			int blocks_to_use = private_key[i];
			if (block1 != 0 && block2 != 0)
			{
				if (i % 2 == 0) //even, then first vector is used
				{
					int item_per_block = ceil(float(first.size()) / float(tmp_blk_one));
					if (blocks_to_use <= block1) // blocks to use < available blocks
					{
						int total_items = item_per_block * blocks_to_use; //total elements
						total_items = (total_items < first_ele) ? (total_items) : (first_ele);
						//int item_count_check = 0; //count the items continuously
						for (int k = 0; k < total_items; k++)
						{
							if (cell_fill % matrix_dimen == 0) //means a column has been filled completely
							{
								Matrix[row_to_fill].push_back(*pointFirst++);
								row_to_fill++;
							}
							else
							{
								Matrix[row_to_fill].push_back(*pointFirst++);
							}
							cell_fill++;
						}
						block1 = block1 - blocks_to_use; //used some blocks
					}
					else
					{
						int total_items = item_per_block * block1; //total elements
						total_items = (total_items < first_ele) ? (total_items) : (first_ele);
						//int item_count_check = 0; //count the items continuously
						for (int k = 0; k < total_items; k++)
						{
								if (cell_fill % matrix_dimen == 0) //means a column has been filled completely
								{
									Matrix[row_to_fill].push_back(*pointFirst++);
									row_to_fill++;
								}
								else
								{
									Matrix[row_to_fill].push_back(*pointFirst++);
								}
								cell_fill++;
						}
						block1 = 0; //used all the blocks
					}
				}
				else if (i % 2 == 1) //even, then first vector is used
				{
					int item_per_block = ceil(float(second.size()) / float(tmp_blk_two));
					if (blocks_to_use <= block2) // blocks to use < available blocks
					{
						int total_items = item_per_block * blocks_to_use; //total elements
						total_items = (total_items < second_ele) ? (total_items) : (second_ele);
						//int item_count_check = 0; //count the items continuously
						for (int k = 0; k < total_items; k++)
						{
							if (cell_fill % matrix_dimen == 0) //means a column has been filled completely
							{
								Matrix[row_to_fill].push_back(*pointSecond++);
								row_to_fill++;
							}
							else
							{
								Matrix[row_to_fill].push_back(*pointFirst++);
							}
							cell_fill++;
						}
						block2 = block2 - blocks_to_use; //used some blocks
					}
					else
					{
						int total_items = item_per_block * block2; //total elements
						total_items = (total_items < second_ele) ? (total_items) : (second_ele);
						//int item_count_check = 0; //count the items continuously
						for (int k = 0; k < total_items; k++)
						{
							if (cell_fill % matrix_dimen == 0) //means a column has been filled completely
							{
								Matrix[row_to_fill].push_back(*pointSecond++);
								row_to_fill++;
							}
							else
							{
								Matrix[row_to_fill].push_back(*pointSecond++);
							}
							cell_fill++;

						}
						block2 = 0; //used all the blocks
					}

				}
			}
			else if (block1 != 0)
			{
				int item_per_block = ceil(float(first.size()) / float(tmp_blk_one));
				int total_items = item_per_block * block1; //total elements
				total_items = (total_items <= first_ele) ? (total_items) : (first_ele);
				//int item_count_check = 0; //count the items continuously
				for (int k = 0; k < total_items; k++)
				{
					if (cell_fill <= total_items)
					{
						if (cell_fill % matrix_dimen == 0) //means a column has been filled completely
						{
							Matrix[row_to_fill].push_back(*pointFirst++);
							row_to_fill++;
						}
						else
						{
							Matrix[row_to_fill].push_back(*pointFirst++);
						}
						cell_fill++;
					}
				}
				block1 = 0; //used all the blocks
				break;
			}
			else if (block2 != 0)
			{
				int item_per_block = ceil(float(second.size()) / float(tmp_blk_two));
				int total_items = item_per_block * block2; //total elements
				total_items = (total_items <= second_ele) ? (total_items) : (second_ele);
				//int item_count_check = 0; //count the items continuously
				for (int k = 0; k < total_items; k++)
				{
					if (cell_fill % matrix_dimen == 0) //means a column has been filled completely
					{
						Matrix[row_to_fill].push_back(*pointSecond++);
                        row_to_fill++;
					}
					else
					{
						Matrix[row_to_fill].push_back(*pointSecond++);
					}
					cell_fill++;
				}
				block2 = 0; //used all the blocks
				break;
			}
			count++;
		}
		if (count == private_key.size())  // rare case, very very rare
		{
			count = 0;
			i = 0;
			std::cout<<"\nreset found\n";
		}
	}

    if(row_to_fill<(matrix_dimen-1))
    {
        int Diff=matrix_dimen-row_to_fill; //no. of remaining row
        for(int i=0;i<Diff;i++)
        {
            addExtra(row_to_fill,matrix_dimen); //matrix_dimen sent as columns
        }
    }
    else if(row_to_fill==(matrix_dimen-1))
    {
        if(Matrix[row_to_fill].size()!=matrix_dimen) //if columns aren't same
        {
            addExtra(row_to_fill,matrix_dimen);
        }
    }

	return Matrix;
}


vector<int> padd_fill(vector<int> key, int old_size)
{
	int* point = private_key.data(); //points to zero of private key
	int count = 0;
	for (int i = old_size; i < key.size(); i++)
	{
		if (count < private_key.size())
		{
			//puts value in key[i] and increments automatically
			//module is needed in order to get index range within size of key
			//key[i] = key[(*point++) % (key.size())];
			key[i] = key[*point++];
			count++;
		}
		if (count == private_key.size())
		{
			point = private_key.data(); //reassign pointer to the start
			count = 0;
		}
	}
	return key;
}


vector<int> key_to_vector_generator(char* file_name)
{
	vector<int> key;//biometric vector
	fstream face_key(file_name, ios::in); //opening biometric key consisting file n read mode
	if (!face_key)
		exit(1);
	else
	{
		while (!face_key.eof())
		{
			string s;
			face_key >> s;

			if (s != " " || s != "\n")
			{
				stringstream val(s);
				int value = 0;
				val >> value;
				key.push_back(value); // creating a vector with the bio keys
			}
		}
		face_key.close();
	}
	return key;
}

void add_row_col(char ch)
{
    if(ch=='R')
    {
        int last_row=U.size();
        int dimen=U[0].size(); //columns
        U.resize(U.size()+1); //add an extra row

        int *point=private_key.data();
        int count=0;

        for(int i=0;i<dimen;i++)
        {
            //int index = (*point++)%dimen;
            int index = *point++;
            if(count!=private_key.size())
            {
                if(i%2==0)
                    U[last_row].push_back(U[0][index]);
                else
                    U[last_row].push_back(U[index][0]);

                count++;
            }
            else
            {
                count=0;
                point=private_key.data();
                continue;
            }
        }

    }
    else
    {
        int *point=private_key.data();
        int count=0;
        int col_size=U[0].size();
        for(int i=0;i<U.size();i++)
        {
            //int index=(*point++)%col_size;
            int index=*point++;
            count++;
            if(count<private_key.size())
            {
                if(i%2==0)
                    U[i].push_back(U[0][index]);
                else
                    U[i].push_back(U[index][0]);
            }
            else
            {
                count=0;
                point=private_key.data();
                continue;
            }
        }

    }
}

int main()
{
	char* bio_name = "face_key.txt";
	char* pri_name = "prime_key.txt";
	char* private_name = "private_key.txt";

	vector<int> bio_key = key_to_vector_generator(bio_name);//biometric vector
	vector<int> pri_key = key_to_vector_generator(pri_name);//prime key vector
	private_key = key_to_vector_generator(private_name);//private key vector

	long int m = bio_key.size();
	long int n = pri_key.size();
	long int s = m + n; // s = m + n contains both their size
	long int q = ceil(sqrt(s));


	int nz1 = q - (m % q); // padding for bio_key vector
	int nz2 = q - (n % q); // padding for pri_key vector

	bio_key.resize(m + nz1, -1); //resize bio vector to accomodate padding, initially -1, as pixels and key can't have -ve values
	pri_key.resize(n + nz2, -1); //resize prime vector to accomodate padding

	bio_key = padd_fill(bio_key, m); //padding fill
	pri_key = padd_fill(pri_key, n);


	nblock_bio = bio_key.size() / q; //dividing vector into blocks
	nblock_pri = pri_key.size() / q;

	int Pad = nz1 + nz2;
	int old_bio_size = bio_key.size();
	int old_pri_size = pri_key.size();

	bio_key.resize(bio_key.size() + Pad, -1); // adding the padded size to the two vectors
	pri_key.resize(pri_key.size() + Pad, -1);

	bio_key = padd_fill(bio_key, old_bio_size); //padding fill
	pri_key = padd_fill(pri_key, old_pri_size);

	// now which vector will start the matrix 'U' is decided by first component of private key module 2
	int index = (private_key[0]) % 2;
	if (index == 0)
		U = matrix_generation(bio_key.size(), pri_key.size(), bio_key, pri_key, nblock_bio, nblock_pri); // even
	else
		U = matrix_generation(pri_key.size(), bio_key.size(), pri_key, bio_key, nblock_pri, nblock_bio); // odd



	// this here makes the matrix dimension divisible by 3, till here we will obtain a sq matrix (probably)
	int row=U.size();
	int col=U[0].size();
    if(row%3!=0) //if not divisible by 3, then add remaining rows
	{
		int time=3-(row%3);
		for(int i=0;i<time;i++)
            add_row_col('R');
	}
    if(col%3!=0)
	{
		int time=3-(col%3);
		for(int i=0;i<time;i++)
            add_row_col('C');
	}
    permutation_and_multiplication(); //for permutation and post-multiplication of U and permutedMatrix


    hybrid_generation();



         for(int j=0;j<private_key.size();j++)
            cout<<private_key[j]<<" ";
        cout<<endl<<endl;
        for(int j=0;j<bio_key.size();j++)
            cout<<bio_key[j]<<" ";

    /*
        for(int i=0;i<HybridMatrix.size();i++)
    {
         for(int j=0;j<HybridMatrix[i].size();j++)
            cout<<HybridMatrix[i][j]<<" ";
        cout<<"\n";
    }
    */

	return 0;
}

