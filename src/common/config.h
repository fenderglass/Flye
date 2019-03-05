//(c) 2016 by Authors
//This file is a part of ABruijn program.
//Released under the BSD license (see LICENSE file)

#pragma once

#include <unistd.h>
#include "flat_hash_map/flat_hash_map.hpp"
#include <fstream>
#include "logger.h"
#include "utils.h"

namespace
{
	std::string trimString(const std::string& str)
	{
		size_t left = 0;
		size_t right = str.size() - 1;
		while(left < str.size() && std::isspace(str[left])) ++left;
		while(right > 0 && std::isspace(str[right])) --right;

		if (right + 1 > left) return str.substr(left, right - left + 1);
		return str;
	}
}

class Config
{
public:
	static Config& instance()
	{
		static Config cfg;
		return cfg;
	}

	static void load(const std::string& filename)
	{
		std::ifstream fin(filename);
		if (!fin) throw std::runtime_error("Can't open config file: " + filename);

		std::string buffer;
		Logger::get().debug() << "Parameters:";
		while(!fin.eof())
		{
			std::getline(fin, buffer);
			if (buffer.empty() || buffer.front() == '#') continue;

			auto tokens = splitString(buffer, '=');
			if (tokens.size() != 2) 
			{
				throw std::runtime_error("Error parsing config file");
			}
			std::string key = trimString(tokens[0]);
			std::string value = trimString(tokens[1]);
			Config::instance()._parameters[key] = std::atof(value.c_str());
			Logger::get().debug() << "\t" << key << "=" << value;
		}
	}

	static float get(const std::string& key) 
	{
		auto itVal = Config::instance()._parameters.find(key);
		if (itVal == Config::instance()._parameters.end())
		{
			throw std::runtime_error("No such parameter: " + key);
		}
		return itVal->second;
	}

private:
	Config(){}
	ska::flat_hash_map<std::string, float> _parameters;
};

struct Parameters
{
	static Parameters& get()
	{
		static Parameters param;
		return param;
	}

	int 	minimumOverlap;
	size_t 	kmerSize;
	size_t 	numThreads;
	bool 	unevenCoverage;
};
