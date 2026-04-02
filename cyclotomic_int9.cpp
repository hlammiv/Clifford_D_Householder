#include "cyclotomic_int9.h"

	/* *********************** */
	/* CONSTRUCTORS  */
	/* *********************** */
	ringZ9::ringZ9(int arr[9]){
		
		for(int i = 0; i < 9; ++i) {
			element[i] = arr[i];
		}
		
		reduce();
	}
	

	ringZ9::ringZ9(int n){
		element[0] = n;
		
		for(int i = 1; i < 9; i++){
			element[i] = 0;
		}
	}
	

	ringZ9::ringZ9(int coeff, int term){
		
		for(int i = 0; i < 9; i++){
			if(i == term%9){
				element[i] =coeff;
			} else {
				element[i] = 0;
			}
		}
		
		reduce();
	}
	
	ringZ9::ringZ9(){

		for(int i = 0; i < 9; i++){
				element[i] = 0;
		}
	}
	
	/* ******* */
	/* INFO */
	/* ******* */

	int ringZ9::getTerm(int i) const {
		if(i > 8 || i < 0){
			return 0;
		} else {
			return element[i];
		}	
	}
	
	

	array<int, 6> ringZ9::getStdArray() const{
		
		array<int, 6> arr;
		
		for(int i = 0; i < 6; i++){
			arr[i] = element[i];
		}
		return arr;
	}


	/* *************** */
	/* MUTATORS */
	/* *************** */

	void ringZ9::reduce(){

		if(isZero()){
			for(int i = 0; i < 9; i++){
				element[i] = 0;
			}
			return;
		}

		for(int i = 6; i < 9; i++){
			for(int j = 1; j < 3; j++){
				element[(i + 3*j)%9] = element[(i + 3*j)%9] - element[i];
			}
		}
		element[6] = 0;
		element[7] = 0;
		element[8] = 0;
	}
	
	
	void ringZ9::scalar_mult(int scalar){
		for(int i = 0; i < 9; i++){
			element[i] = scalar*element[i];
		}
	}
	
	void ringZ9::scalar_div(int scalar){
		if(scalar == 0){
			cout << "Zero passed as scalar in ringZ9::scalar_div(int scalar)" << endl;
		}else {		
			for(int i = 0; i < 9; i++){
				element[i] = element[i]/scalar;
			}
		}
	}
	
	
	
	/* ***************** */
	/* OPERATIONS */
	/* ***************** */
	ringZ9 ringZ9::operator+(const ringZ9& right) const{
		ringZ9 sum;
		
		for(int i = 0; i < 6; i++){
			sum.element[i] = element[i] + right.element[i];
		}

		return sum;
	}
	
	
	ringZ9 ringZ9::operator-(const ringZ9& right) const{
		ringZ9 diff;
		
		for(int i = 0; i < 6; i++){
			diff.element[i] = element[i] - right.element[i];
		}

		return diff;
	}
	
	ringZ9 ringZ9::operator*(const ringZ9& right) const{
		ringZ9 prod;
		
		for(int i = 0; i <6; i++){
			for(int j = 0; j < 6; j++){
				prod.element[(i + j)%9] = prod.element[(i + j)%9] + element[i]*right.element[j];
			}
		}
		
		prod.reduce();
		return prod;
	}
	

	ringZ9 ringZ9::operator*(const int& right_scalar) const{
		ringZ9 prod;
		
		for(int i = 0; i < 6; i++){
			prod.element[i] = element[i]*right_scalar;
		}
		
		return prod;
	}
	

	ringZ9 ringZ9::operator/(const int& right_scalar) const{
		ringZ9 prod;
		
		if(right_scalar == 0){
			cout << "Zero passed as scalar in ringZ9::operator/(const int& right_scalar)" << endl;
			return prod;
		}
		
		for(int i = 0; i < 6; i++){
			prod.element[i] = element[i]/right_scalar;
		}
		
		return prod;
	}
	
	
	ringZ9 ringZ9::operator=(const ringZ9& right) {
		if (this != &right) {
			for (int i = 0; i < 9; ++i)
				this->element[i] = right.element[i];
			}
		return *this;
	}
	
	
	bool ringZ9::operator==(const ringZ9& right) const{
		// Direct coefficient comparison: avoids allocating a temporary ringZ9
		// and performing a full subtraction just to call isZero().
		// Both sides are already in reduced form (elements[6..8] == 0), so
		// comparing the first 6 coefficients is sufficient.
		for (int i = 0; i < 6; ++i)
			if (element[i] != right.element[i]) return false;
		return true;
	}
	
	bool ringZ9::operator!=(const ringZ9& right) const{
		for (int i = 0; i < 6; ++i)
			if (element[i] != right.element[i]) return true;
		return false;
	}
	
	
	/* *************************** */
	/* TYPE CONVERSIONS */
	/* *************************** */
	double ringZ9::real_part() const {
                static const double cos_vals[4] = {
                0.766044443118978, 		//cos(2 * M_PI / 9)
                0.17364817766693041, 	//sin(M_PI / 18),
                -0.5, 					//cos(6 pi/9)
                -0.9396926207859083, 	//-cos(M_PI / 9),
                };
                double sum = element[0];
                sum += cos_vals[0] * element[1];
                sum += cos_vals[1] * element[2];
                sum += cos_vals[2] * element[3];
                sum += cos_vals[3] * (element[4]+element[5]);
                return sum;
        }
        
        double ringZ9::imag_part() const {
                static const double sin_vals[4] = {
                        0.6427876096865393,  	// sin(2π/9)
                        0.984807753012208,   	// sin(4π/9)
                        0.8660254037844387,   	// sin(6π/9)
                        0.3420201433256689,   	// sin(8π/9)
                };
                double sum = 0.0;
                sum += sin_vals[0] * element[1];
                sum += sin_vals[1] * element[2];
                sum += sin_vals[2] * element[3];
                sum += sin_vals[3] * (element[4]-element[5]); //The minus sign is due to periodicity
                return sum;
        }
	
	double ringZ9::complexArg() const{
		return atan2(imag_part(), real_part());
	}
	
	
	double ringZ9::abs_val() const{
		return sqrt(real_part()*real_part() + imag_part()*imag_part());
	}
	
	
	double ringZ9::abs_val_sq() const{
		return real_part()*real_part() + imag_part()*imag_part();
	}
	
	
	complex<double> ringZ9::toComplexDouble() const{
		complex<double> z(real_part(),imag_part());
		return z;
	}
	
	
	/* **************************** */
	/* NUMBER THEORETIC */
	/* **************************** */
	// Complex conjugate
	ringZ9 ringZ9::complexConj() const{
		int swap[9] = {0,0,0,0,0,0,0,0,0};
		
		for (int i = 0; i < 6; i++){
			swap[(8*i)%9] = element[i];
		}
		return ringZ9(swap);
	}
	
	

	ringZ9 ringZ9::GaloisAut(int k) const{
		k = k%9;
		int swap[9] = {0,0,0,0,0,0,0,0,0};
		
		for (int i = 0; i < 6; i++){
			swap[(((k*i % 9) + 9) % 9)] = swap[(((k*i % 9) + 9) % 9)] + element[i];
		}
		return ringZ9(swap);
	}
	

	int ringZ9::fieldNorm() const{
		ringZ9 prod(1);
		
		for(int i = 1; i < 9; i++){
			if(gcd(i,9) == 1){
				prod = prod*GaloisAut(i);
			}
		}
		
		prod.reduce();
 		return prod.getTerm(0);
	}
	

	int ringZ9::fieldTrace() const{
		ringZ9 sum = GaloisAut(1);
		
		for(int i = 2; i < 9; i++){
			if(gcd(i,9) == 1){
				sum = sum+GaloisAut(i);
			}
		}
 		return sum.getTerm(0);
	}
	
	
	ringZ9 ringZ9::partialFieldNorm() const{
		ringZ9 prod = GaloisAut(2);
		
		for(int i = 4; i < 9; i++){
			if(gcd(i,9) == 1){
				prod = prod*GaloisAut(i);
			}
		}
		return prod;		
	}
	
	int ringZ9::tauFieldNorm() const{
		ringZ9 prod = GaloisAut(1)*GaloisAut(2)*GaloisAut(5);
		
		return prod.getTerm(0);
	}
	

	int ringZ9::weight() const{
		int sum = 0;
		
		for(int i = 0; i < 6; i++){
			sum = sum + abs(element[i]);
		}
		return sum;
	}
	
	int ringZ9::signedWeight() const{
		int sum = 0;
		
		for(int i = 0; i < 6; i++){
			sum = sum + element[i];
		}
		return sum;
	}
	
	int ringZ9::quad() const{
		int sum = 0;
		
		for(int i = 0; i < 6; i++){
			sum = sum + element[i]*element[i];
		}
		
		return sum - element[0]*element[3] - element[1]*element[4] - element[2]*element[5];
	}
	
	
	ringZ9 ringZ9::formalDerivative() const{
		
		int new_element[9] = {0};
		
		for(int i = 0; i < 6; i++){
			new_element[i] = (i+1)*element[i+1];
		}
		
		return ringZ9(new_element);
	}
	
	
	int ringZ9::sdeChi() const{
		if(signedWeight() % 3 != 0) return 0;
		
		ringZ9 test = formalDerivative();
		
		for(int i = 1; i < 6; i++){
			if((test.signedWeight() / i)% 3 != 0 ) return i;
			
			test = (test.formalDerivative())/i;
			
		}
		
		
		return 6;
	}
	
	
	/* ********* */
	/* TESTS */
	/* ********* */
	bool ringZ9::isZero() const{
		
		for(int i = 0; i < 3; i++){
			for(int j = 1; j < 3; j++){
				
				if(element[i + 3*j] !=element[i]){
					return false;
				}
			}
		}
		
		return true;
	}
	
	
	bool ringZ9::isReal() const{
		ringZ9 test = *this - this->complexConj();
		return test.isZero();
	}
	
	
	bool ringZ9::isImag() const{
		ringZ9 test = *this + this->complexConj();
		return test.isZero();
	}
	

	bool ringZ9::isInt() const{
		for(int i = 1; i < 6; i++){
			if(element[i] != 0){
				return false;
			}
		}
		return true;
	}
	
	bool ringZ9::isDivisibleByInt(int k) const{
		for(int i = 0; i < 6; i++){
			if(element[i]%k != 0){
				return false;
			}
		}
		return true;
	}
	
	
	/* ********** */
	/* PRINTS */
	/* ********** */
	void ringZ9::print() const{
		
		if(isZero()){ // prints zero if element equals zero
			cout << "(0)" << endl;
			return;
		}
		
		int i = 0;
		int j = 8;
		
		while(element[i] == 0){ i++;} // i is now the exponent on the first nonzero term
		while(element[j] == 0){	j--;} // j is now the expoent on the last nonzero term
		
		
		if(i == j){
			cout << "(" << element[i] << "z^" << i << ")";
			return;
		}
		
		cout << "(";
		
		for(int k = i; k < j; k++){ 
			if(element[k] != 0){
				cout << element[k] << "z^" << k << " + ";
			}
		}
		
		cout << element[j] << "z^" << j << ")";
	}
	
