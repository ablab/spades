////***************************************************************************
////* Copyright (c) 2011-2014 Saint-Petersburg Academic University
////* All Rights Reserved
////* See file LICENSE for details.
////****************************************************************************
// todo remove!!!
//#ifndef IS_CORRUPTING_WRAPPER_HPP_
//#define IS_CORRUPTING_WRAPPER_HPP_
//
//namespace io {
//
//class ISCorruptingWrapper: public DelegatingReaderWrapper<PairedRead> {
//private:
//	const size_t is_;
//public:
//	typedef PairedRead ReadType;
//
//	explicit ISCorruptingWrapper(IReader<ReadType>& reader, size_t is) :
//			DelegatingReaderWrapper<PairedRead>(reader), is_(is) {
//	}
//
//	/* virtual */
//	ISCorruptingWrapper& operator>>(ReadType& read) {
//		(this->reader()) >> read;
//		read = PairedRead(read.first(), read.second(), is_);
//		return *this;
//	}
//
//};
//
//}
//
//#endif /* IS_CORRUPTING_WRAPPER_HPP_ */
