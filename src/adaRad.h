#define RADEMAX 15
int bucket[(1 << RADEMAX)];
static inline  int adaRad(int n)
{
	// const  int cute[]={2,3,4,5,6,7,8,9,10,11,12,13,14};
	const  int cutn[]={19,75,300,1208,4860,19551,78657,316442,1273073,5121682,20604957,82895478,333495489};
	if (n > 333495489) return (int)(RADEMAX);
	int * p = (int *)(cutn); 
	while (n>*(p++));
	return (int)(p-cutn+1);
}
