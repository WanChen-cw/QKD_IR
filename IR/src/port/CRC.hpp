#ifndef CRC_HPP_
#define CRC_HPP_

#include <vector>


template <typename B = int>
class CRC
{
private:
	uint32_t crcPolynomial; //0x04C11DB7 0x814141AB 0x04C11DB7 320x1EDC6F41 0x814141AB 0x32583499
public:
	CRC(const uint32_t crcPolynomial) :crcPolynomial(crcPolynomial){}

	~CRC() = default;
	
	uint32_t calculate(const std::vector<B>& data);

	std::vector<uint32_t> segcrccalculate(const std::vector<B>& U, std::vector<int>& SegSite);

	uint32_t sitecrccalculate(const std::vector<B>& U,int site);
};

template<typename B>
uint32_t CRC<B>::calculate(const std::vector<B>& data){
	uint32_t crc = 0xFFFFFFFF;
	for (auto byte : data) {
		crc ^= (byte << 24);
		for (int j = 0; j < 8; ++j) {
			if (crc & 0x80000000)
				crc = (crc << 1) ^ crcPolynomial;
			else
				crc <<= 1;
		}
	}
	return crc;
}

template<typename B>
std::vector<uint32_t> CRC<B>
::segcrccalculate(const std::vector<B>& U, std::vector<int> &SegSite){
	std::vector<uint32_t> crc_seg;
	uint32_t crc = 0xFFFFFFFF;
	int idx_segsite = 0;
	for (size_t i = 0; i < U.size(); ++i) {
	    crc ^= (U[i] << 24);
	    for (int j = 0; j < 8; ++j) {
	        if (crc & 0x80000000)
	            crc = (crc << 1) ^ crcPolynomial;
	        else
	            crc <<= 1;
	    }
		if (idx_segsite < SegSite.size()) {
			if (i == SegSite[idx_segsite]) {
				crc_seg.push_back(crc);
				idx_segsite++;
			}
		}
	}
	return crc_seg;
}

template<typename B>
uint32_t CRC<B>
::sitecrccalculate(const std::vector<B>& U, int site){
	uint32_t crc = 0xFFFFFFFF;
	for (size_t i = 0; i <= site; ++i) {
		crc ^= (U[i] << 24);
		for (int j = 0; j < 8; ++j) {
			if (crc & 0x80000000)
				crc = (crc << 1) ^ crcPolynomial;
			else
				crc <<= 1;
		}
	}
	return crc;
}

#endif /*CRC_HPP_ */


