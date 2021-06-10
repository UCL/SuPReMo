// ====================================================================================================   
//                                                                                                        
//   SuPReMo: Surrogate Parameterised Respiratory Motion Model                                            
//            An implementation of the generalised motion modelling and image registration framework      
//                                                                                                        
//   Copyright (c) University College London (UCL). All rights reserved.                                  
//                                                                                                        
//   This software is distributed WITHOUT ANY WARRANTY; without even                                      
//   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR                                  
//   PURPOSE.                                                                                             
//                                                                                                        
//   See LICENSE.txt in the top level directory for details.                                              
//                                                                                                        
// ====================================================================================================   





//----------
// Includes
//----------
#include "CommandLineParser.h"
#include <algorithm>
#include <iomanip> 

//------------------------
// splitStringbyDelimiter
//------------------------
std::vector<std::string> splitStringbyDelimiter(const std::string & stringToSplitIn, const std::string & delimiter)
{
	// Copy the input since we need to erase bist of this string 
	std::string stringToSplit(stringToSplitIn);

	// The output string
	std::vector<std::string> vSplitString;

	// Position 
	size_t pos = 0;
	
	while ((pos = stringToSplit.find(delimiter)) != std::string::npos) {
		if (!stringToSplit.substr(0, pos).empty())
		{
			vSplitString.push_back(stringToSplit.substr(0, pos));
		}
		stringToSplit.erase(0, pos + delimiter.length());
	}

	// Only append item if it is not empty. 
	if (!stringToSplit.empty())
	{
		vSplitString.push_back(stringToSplit);
	}

	return vSplitString;
}


//----------------------------------
// printFormattedCommandLineOptions
//----------------------------------
void printFormattedCommandLineOptions(const std::map<std::string, CommandLineOption>& commandLineOptions,
  const std::string& optionSectionHeading, const unsigned int maxWidth, const std::vector<std::string>& optionsToPrint) 
{
  // Print the heading of this section
  std::cout << optionSectionHeading << std::endl << std::string(optionSectionHeading.size(), '~') << std::endl;
  
  std::vector<std::string> optionsToPrintAndInCommandLineOptions;

  // Make sure all required options that should be printed are present in the command line options
  for (auto optionIt = optionsToPrint.begin(); optionIt != optionsToPrint.end(); ++optionIt)
  {
    if (commandLineOptions.find(*optionIt) == commandLineOptions.end())
      std::cout << "WARNING: Option required to print: " << *optionIt << " but not within command line options!" << std::endl;
    else 
      optionsToPrintAndInCommandLineOptions.push_back(*optionIt);
  }
  
  // Determine the maximum length of the option and format combined
  size_t optionAndFormatWidth = 0;
  for (auto optionIt = optionsToPrintAndInCommandLineOptions.begin(); optionIt != optionsToPrintAndInCommandLineOptions.end(); ++optionIt)
  {
    commandLineOptions.find(*optionIt)->second.format.size();
    optionAndFormatWidth = (optionIt->size() + commandLineOptions.find(*optionIt)->second.format.size()) > optionAndFormatWidth ?
      optionIt->size() + commandLineOptions.find(*optionIt)->second.format.size() : optionAndFormatWidth;
  }
  optionAndFormatWidth += 3;
  
  if (maxWidth < optionAndFormatWidth)
  {
    std::cout << "Cannot print parameters, increase maximum print width." << std::endl;
    return;
  }

  size_t descriptionWidth = maxWidth - optionAndFormatWidth;
  
  for (auto optionIt = optionsToPrintAndInCommandLineOptions.begin(); optionIt != optionsToPrintAndInCommandLineOptions.end(); ++optionIt)
  {
    // Split the description into paragraphs
    auto paragraphs = splitStringbyDelimiter(commandLineOptions.find(*optionIt)->second.description, "\n");
    
    // Track if we are printing the first help line
    bool firstLine = true;

    for (size_t i = 0; i < paragraphs.size(); ++i)
    {
      // split paragraph into words
      auto wordsInParagraph = splitStringbyDelimiter(paragraphs[i], " ");
      
      // Catch the case where there was no description provided...
      if (wordsInParagraph.empty())
      {
        wordsInParagraph.push_back(" ");
      }

      while (!wordsInParagraph.empty())
      {
        // Add first word regardless even if it goes beyond description length...
        std::string curParagraphLine = wordsInParagraph[0];
        wordsInParagraph.erase(wordsInParagraph.begin());

        // Add words to the current paragraph line until the description width was reached
        while ((!wordsInParagraph.empty() && (curParagraphLine.size() + 1 + wordsInParagraph[0].size() < descriptionWidth)))
        {
          curParagraphLine += " " + wordsInParagraph[0];
          wordsInParagraph.erase(wordsInParagraph.begin());
        }

        if (firstLine)
        {
          firstLine = false;
          std::string s = *optionIt + " " + commandLineOptions.find(*optionIt)->second.format;;
          std::cout << std::setw(optionAndFormatWidth) << std::left << s << curParagraphLine << std::endl;
        }
        else
        {
          std::cout << std::setw(optionAndFormatWidth) << std::left << " " << curParagraphLine << std::endl;
        }
      }
    }
  }
  std::cout << std::endl;
  return;
}



//--------------------------------------
// CommandLineParser::CommandLineParser
//--------------------------------------
CommandLineParser::CommandLineParser(int &argc, char **argv, std::map<std::string, CommandLineOption> & allowedCommandLineOptionsIn)
    :allRequiredParametersSet(true)
{
  // Save the executable name
  this->executableName = std::string(argv[0]);

  // And every other item from the command line
  for ( int i = 1; i < argc; ++i )
  {
    this->tokens.push_back( std::string( argv[i] ) );
  }
    
  // Check that all tokens are valid as defined by the allowed command line options
  // that is 
  // 1) check that the option is known 
  // 2) skip the N further components of that option before proceeding
  bool* knownAndChecked = new bool[this->tokens.size()]; // ensure initialisation to false
  for ( int i = 0; i < this->tokens.size(); ++i )
    knownAndChecked[i] = false;

  for ( int i = 0; i < this->tokens.size(); ++i )
  {
    // find the token in the map of known parameters
    auto curCMDOption = allowedCommandLineOptionsIn.find( this->tokens[i] );
      
    if ( curCMDOption == allowedCommandLineOptionsIn.end() )
    {
      char msg[200];
      sprintf_s( msg, "Command line option: %s NOT KNOWN", this->tokens[i].c_str() );
      supremo_print_error( msg );
      supremo_exit( 1, __FILE__, __LINE__ );
    }
    else{
      // set the known and checked for the parameter to true and also the n numbers required by that command
      knownAndChecked[i] = true;
      int curNumOfParameterParts = allowedCommandLineOptionsIn[this->tokens[i]].numberOfComponents;
      for ( int j = 0; j < curNumOfParameterParts; ++j )
      {
        ++i;
        knownAndChecked[i] = true;
      }
    }
  }

  // Make sure that all required parameters were provided
  typedef std::map<std::string, CommandLineOption>::const_iterator AllowedCommandLineOptsIteratorType;
  for ( AllowedCommandLineOptsIteratorType iter = allowedCommandLineOptionsIn.begin(); 
        iter != allowedCommandLineOptionsIn.end(); 
        ++iter )
  {
    if ( iter->second.required && (! this->cmdOptionExists(iter->first)) )
    {
      char msg[200];
      sprintf_s( msg, "Required command line option: %s was not provided", iter->first.c_str() );
      supremo_print_warning( msg );
      this->allRequiredParametersSet = false;
    }
  }
}




//-----------------------------------------
// CommandLineParser::getCmdOptionAsString
//-----------------------------------------
const std::string& CommandLineParser::getCmdOptionAsString( const std::string &option, int skip ) const
{
  std::vector<std::string>::const_iterator itr;
  itr = std::find(this->tokens.begin(), this->tokens.end(), option);

  // if option was not found return empty string
  if (itr == tokens.end())
  {
    static const std::string empty_string("");
    return empty_string;
  }

  // skip forward if requested
  for (int i = 0; i < skip; ++i)
  {
    ++itr;
    if (itr == tokens.end())
    {
      static const std::string empty_string("");
      return empty_string;
    }
  }

  // return the requested option if possible
  if (itr != this->tokens.end() && ++itr != this->tokens.end()){
    return *itr;
  }
    
  // otherwise default to empty string
  static const std::string empty_string("");
  return empty_string;
}




//--------------------------------------
// CommandLineParser::getCmdOptionAsInt
//--------------------------------------
const int CommandLineParser::getCmdOptionAsInt( const std::string &option, int skip ) const
{
  std::string optionString = this->getCmdOptionAsString( option, skip );
  if (!optionString.empty())
  {
    return atoi( optionString.c_str() );
  }
  else
  {
    char msg[200];
    sprintf_s( msg, "Command line option: %s COULD NOT BE CONVERTED TO INT.", optionString.c_str() );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }
  return 0;
}




//----------------------------------------
// CommandLineParser::getCmdOptionAsFloat
//----------------------------------------
const float CommandLineParser::getCmdOptionAsFloat( const std::string &option, int skip ) const
{
  std::string optionString = this->getCmdOptionAsString( option, skip );
  if (!optionString.empty())
  {
    return static_cast<float>(atof( optionString.c_str() ));
  }
  else
  {
    char msg[200];
    sprintf_s( msg, "Command line option: %s COULD NOT BE CONVERTED TO FLOAT.", optionString.c_str() );
    supremo_print_error( msg );
    supremo_exit( 1, __FILE__, __LINE__ );
  }
  return 0;
}




//------------------------------------
// CommandLineParser::cmdOptionExists
//------------------------------------
bool CommandLineParser::cmdOptionExists( const std::string &option ) const
{
  return std::find(this->tokens.begin(), this->tokens.end(), option)
    != this->tokens.end();
}
  



//-----------------------------------
// CommandLineParser::getCommandLine
//-----------------------------------
std::string CommandLineParser::getCommandLine() const
{
  // Append the executable name
  std::string completeCommandLine = executableName;

  // Add each command line token
  for (size_t i = 0; i < this->tokens.size(); ++i)
  {
    completeCommandLine.append(" ");
    completeCommandLine.append(this->tokens[i]);
  }

  return completeCommandLine;      
}
