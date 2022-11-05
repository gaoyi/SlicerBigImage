/*
 * utils.h
 *
 *  Created on: Jul 15, 2011
 *      Author: tcpan
 */
#include "Logger.h"
#include <sstream>
#include <fstream>
#include <iostream>
#include <string.h>
#include <cstdlib>
#include <iomanip>

namespace cci {
namespace common {

const int event::COMPUTE = 0;
const int event::MEM_IO = 11;
const int event::GPU_MEM_IO = 12;
const int event::NETWORK_IO = 21;
const int event::NETWORK_WAIT = 22;
const int event::NETWORK_IO_NB = 23;
const int event::FILE_I = 31;
const int event::FILE_O = 32;
const int event::ADIOS_INIT = 41;
const int event::ADIOS_OPEN = 42;
const int event::ADIOS_ALLOC = 43;
const int event::ADIOS_WRITE = 44;
const int event::ADIOS_CLOSE = 45;
const int event::ADIOS_BENCH_OPEN = 52;
const int event::ADIOS_BENCH_ALLOC = 53;
const int event::ADIOS_BENCH_WRITE = 54;
const int event::ADIOS_BENCH_CLOSE = 55;
const int event::ADIOS_FINALIZE = 46;
const int event::OTHER = -1;

std::string event::getAsString() const {
	std::stringstream ss;
	ss << std::fixed << "[" << id << "]" << eventtype << "=\t" << name << ":\t" << starttime;
	if (endtime != -1)
		ss << "\t-\t" << endtime << "\t=\t" << (endtime - starttime + 1);
	else
		ss << "\t\t\t\t";
	ss << "\t" << annotation;
	return ss.str();
}


LogSession::LogSession(const int &_id, const std::string &_name, const int &_group,
		const std::string &_session_name, long long &_start) :
	id(_id), name(_name), group(_group), session_name(_session_name), events(), start(_start) {
	events.clear();

	countByEventName.clear();
	sumDurationByEventName.clear();
	sumSquareDurationByEventName.clear();;

	countByEventType.clear();
	sumDurationByEventType.clear();
	sumSquareDurationByEventType.clear();;
};
LogSession::~LogSession() {
	events.clear();
	countByEventName.clear();
	sumDurationByEventName.clear();
	sumSquareDurationByEventName.clear();;

	countByEventType.clear();
	sumDurationByEventType.clear();
	sumSquareDurationByEventType.clear();;
};
void LogSession::restart() {
	events.clear();
	countByEventName.clear();
	sumDurationByEventName.clear();
	sumSquareDurationByEventName.clear();;

	countByEventType.clear();
	sumDurationByEventType.clear();
	sumSquareDurationByEventType.clear();;
}
void LogSession::log(cci::common::event e) {
	events.push_back(e);

	std::string ename = e.getName();
	int etype = e.getType();
	long long duration = e.getEnd() - e.getStart() + 1;
	countByEventName[ename] += 1;
	sumDurationByEventName[ename] += duration;
	sumSquareDurationByEventName[ename] += duration * duration;
	countByEventType[etype] += 1;
	sumDurationByEventType[etype] += duration;
	sumSquareDurationByEventType[etype] += duration * duration;
};

void LogSession::toString(std::string &header, std::string &value) {
	std::stringstream ss1, ss2;
	ss1 << "pid,hostName,group,sessionName,";
	ss2 << id << "," << name << "," << group << "," << session_name << "," << std::fixed;

	for (int i = 0; i < events.size(); ++i) {
		ss1 << events[i].getName() << "," << events[i].getType() << ",";
		ss2 << (events[i].getStart() - start + 1) << "," << (events[i].getEnd() - start + 1) << ",";
	}
	header.assign(ss1.str());
	value.assign(ss2.str());
};

void LogSession::toOneLineString(std::string &value) {
	std::stringstream ss1;
	ss1 << "pid," << id << ",hostName," << name << ",group," << group << ",sessionName," << session_name << "," << std::fixed;

	for (int i = 0; i < events.size(); ++i) {
		ss1 << events[i].getName() << "," << events[i].getType() << "," << (events[i].getStart() - start + 1) << "," << (events[i].getEnd() - start + 1) << "," << events[i].getAnnotation() << ",";
	}
	value.assign(ss1.str());
};

void LogSession::toSummaryStringByName(std::string &header, std::string &value) {
	std::stringstream ss1, ss2;
	ss1 << "pid,hostName,group,sessionName,";
	ss2 << id << "," << name << "," << group << "," << session_name << "," << std::fixed;

	for (std::tr1::unordered_map<std::string, long>::iterator iter= countByEventName.begin();
			iter != countByEventName.end(); ++iter) {
		std::string name = iter->first;
		ss1 << name << " count," << name << " mean," << name << " stdev, ";
		ss2 << iter->second << "," << getMeanByEventName(name) << "," << getStdevByEventName(name) << ",";
	}
	header.assign(ss1.str());
	value.assign(ss2.str());
}
void LogSession::toSummaryStringByType(std::string &header, std::string &value) {
	std::stringstream ss1, ss2;
	ss1 << "pid,hostName,group,sessionName,";
	ss2 << id << "," << name << "," << group << "," << session_name << "," << std::fixed;

	for (std::tr1::unordered_map<int, long>::iterator iter = countByEventType.begin();
			iter != countByEventType.end(); ++iter) {
		int type = iter->first;
		ss1 << "type " << type << " count," << "type " << type << " mean," << "type " << type << " stdev, ";
		ss2 << iter->second << "," << getMeanByEventType(type) << "," << getStdevByEventType(type) << ",";
	}
	header.assign(ss1.str());
	value.assign(ss2.str());
}


// session id is something like a filename or image name or hostname.
LogSession* Logger::getSession(const std::string &session_name) {
	if (values.find(session_name) == values.end()) {
		cci::common::LogSession session(id, name, group, session_name, starttime);
		values[session_name] = session;
	}
	return &(values[session_name]);
}

std::vector<std::string> Logger::toStrings() {
	// headers
	std::vector<std::string> output;

	for (std::tr1::unordered_map<std::string, cci::common::LogSession >::iterator iter = values.begin();
			iter != values.end(); ++iter) {
		std::string headers;
		std::string times;

		iter->second.toString(headers, times);

		output.push_back(headers);
		output.push_back(times);
	}
	return output;
}
std::vector<std::string> Logger::toOneLineStrings() {
	// headers
	std::vector<std::string> output;

	for (std::tr1::unordered_map<std::string, cci::common::LogSession >::iterator iter = values.begin();
			iter != values.end(); ++iter) {
		std::string times;

		iter->second.toOneLineString(times);

		output.push_back(times);
	}
	return output;
}
std::vector<std::string> Logger::toSummaryStringsByName() {
	// headers
	std::vector<std::string> output;

	for (std::tr1::unordered_map<std::string, cci::common::LogSession >::iterator iter = values.begin();
			iter != values.end(); ++iter) {
		std::string headers;
		std::string times;

		iter->second.toSummaryStringByName(headers, times);

		output.push_back(headers);
		output.push_back(times);
	}
	return output;
}
std::vector<std::string> Logger::toSummaryStringsByType() {
	// headers
	std::vector<std::string> output;

	for (std::tr1::unordered_map<std::string, cci::common::LogSession >::iterator iter = values.begin();
			iter != values.end(); ++iter) {
		std::string headers;
		std::string times;

		iter->second.toSummaryStringByType(headers, times);

		output.push_back(headers);
		output.push_back(times);
	}
	return output;
}
void Logger::write(const std::string &prefix) {
	std::vector<std::string> timings = this->toOneLineStrings();
	std::stringstream ss;
	for (int i = 0; i < timings.size(); i++) {
			ss << timings[i] << std::endl;
	}

	std::stringstream fss;
	fss << prefix << "-" << id << ".csv";

	time_t now;
	time(&now);
	struct tm *current;
	current = localtime(&now);
	struct tm *start;
	start = localtime(&start_t);


	std::ofstream ofs2(fss.str().c_str());
	ofs2 << "v2.1,logger start,";
	ofs2 << std::setw( 4 ) << std::setfill( '0' ) << (start->tm_year + 1900) << "/";
	ofs2 << std::setw( 2 ) << std::setfill( '0' ) << (start->tm_mon + 1) << "/";
	ofs2 << std::setw( 2 ) << std::setfill( '0' ) << start->tm_mday << " ";
	ofs2 << std::setw( 2 ) << std::setfill( '0' ) << start->tm_hour << ":";
	ofs2 << std::setw( 2 ) << std::setfill( '0' ) << start->tm_min << ":";
	ofs2 << std::setw( 2 ) << std::setfill( '0' ) << start->tm_sec;
	ofs2 << ",logger finish,";
	ofs2 << std::setw( 4 ) << std::setfill( '0' ) << (current->tm_year + 1900) << "/";
	ofs2 << std::setw( 2 ) << std::setfill( '0' ) << (current->tm_mon + 1) << "/";
	ofs2 << std::setw( 2 ) << std::setfill( '0' ) << current->tm_mday << " ";
	ofs2 << std::setw( 2 ) << std::setfill( '0' ) << current->tm_hour << ":";
	ofs2 << std::setw( 2 ) << std::setfill( '0' ) << current->tm_min << ":";
	ofs2 << std::setw( 2 ) << std::setfill( '0' ) << current->tm_sec << std::endl;
	ofs2 << ss.str() << std::endl;
	ofs2.close();
}

void Logger::write() {
	write(logprefix);
}


#if defined (WITH_MPI)
void Logger::writeCollectively(const int &rank, const int &manager_rank, MPI_Comm &comm_world) {
	writeCollectively(logprefix, rank, manager_rank, comm_world);
}
void Logger::writeCollectively(const std::string &prefix, const int &rank, const int &manager_rank, MPI_Comm &comm_world) {
	int size;
	MPI_Comm_size(comm_world, &size);

	// now do a collective io for the log
	std::vector<std::string> timings = this->toOneLineStrings();
	std::stringstream ss;
	for (int i = 0; i < timings.size(); i++) {
			ss << timings[i] << std::endl;
	}
	std::string logstr = ss.str();
	int logsize = logstr.size();

	char *sendlog = (char *)malloc(sizeof(char) * logsize + 1);
	memset(sendlog, 0, sizeof(char) * logsize + 1);
	strncpy(sendlog, logstr.c_str(), logsize);
	ss.str(std::string());

	int *recbuf = NULL;

	if (rank == manager_rank)
		recbuf = (int *) malloc(size * sizeof(int));

		// now send the thing to manager
		//      first gather sizes
		MPI_Gather(&logsize, 1, MPI_INT, recbuf, 1, MPI_INT, manager_rank, comm_world);


		//      then gatherv the messages.
		char *logdata = NULL;
		int * displbuf = NULL;


		if (rank == manager_rank) {
		// then perform exclusive prefix sum to get the displacement and the total length
		displbuf = (int *) malloc(size * sizeof(int));
		displbuf[0] = 0;
			for (int i = 1; i < size; i++) {
					displbuf[i] = displbuf[i-1] + recbuf[i-1];
			}
			int logtotalsize = displbuf[size - 1] + recbuf[size - 1];

			logdata = (char*) malloc(logtotalsize * sizeof(char) + 1);
			memset(logdata, 0, logtotalsize * sizeof(char) + 1);

	}
	MPI_Gatherv(sendlog, logsize, MPI_CHAR, logdata, recbuf, displbuf, MPI_CHAR, manager_rank, comm_world);

	free(sendlog);

	if (rank == manager_rank) {
			free(recbuf);
			free(displbuf);

		std::stringstream fss;
		fss << prefix << ".csv";

		time_t now;
		time(&now);
		struct tm *current;
		current = localtime(&now);
		struct tm *start;
		start = localtime(&start_t);


		std::ofstream ofs2(fss.str().c_str());
		ofs2 << "v2.1,logger start,";
		ofs2 << std::setw( 4 ) << std::setfill( '0' ) << (start->tm_year + 1900) << "/";
		ofs2 << std::setw( 2 ) << std::setfill( '0' ) << (start->tm_mon + 1) << "/";
		ofs2 << std::setw( 2 ) << std::setfill( '0' ) << start->tm_mday << " ";
		ofs2 << std::setw( 2 ) << std::setfill( '0' ) << start->tm_hour << ":";
		ofs2 << std::setw( 2 ) << std::setfill( '0' ) << start->tm_min << ":";
		ofs2 << std::setw( 2 ) << std::setfill( '0' ) << start->tm_sec;
		ofs2 << ",logger finish,";
		ofs2 << std::setw( 4 ) << std::setfill( '0' ) << (current->tm_year + 1900) << "/";
		ofs2 << std::setw( 2 ) << std::setfill( '0' ) << (current->tm_mon + 1) << "/";
		ofs2 << std::setw( 2 ) << std::setfill( '0' ) << current->tm_mday << " ";
		ofs2 << std::setw( 2 ) << std::setfill( '0' ) << current->tm_hour << ":";
		ofs2 << std::setw( 2 ) << std::setfill( '0' ) << current->tm_min << ":";
		ofs2 << std::setw( 2 ) << std::setfill( '0' ) << current->tm_sec << std::endl;
   		ofs2 << logdata << std::endl;
		ofs2.close();

		//printf("%s\n", logdata);
		free(logdata);
	}
}
#endif

}}
