#include <pybind11/pybind11.h>
#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
namespace py = pybind11;

#include <nsearch/FASTA/Reader.h>
#include <nsearch/Database.h>
#include <nsearch/Database/HitWriter.h>
#include <nsearch/Database/GlobalSearch.h>
#include <nsearch/Sequence.h>
#include <nsearch/Alphabet/DNA.h>
#include <nsearch/Alphabet/Protein.h>

#include <Search.h>

#include <string>
#include <memory>
#include <stdexcept>

#include <Common.h>
#include "FileFormat.h"
#include "WorkerQueue.h"

template < typename A >
using QueryWithHits     = std::pair< Sequence< A >, HitList< A > >;

template < typename A >
using QueryWithHitsList = std::deque< QueryWithHits< A > >;

template < typename A >
class QueueItemInfo< QueryWithHitsList< A > > {
public:
  static size_t Count( const QueryWithHitsList< A >& list ) {
    return std::accumulate(
      list.begin(), list.end(), 0,
      []( int sum, const QueryWithHits< A >& q ) { return sum + q.second.size(); } );
  }
};

template < typename A >
class SearchResultsWriterWorker {
public:
  SearchResultsWriterWorker( const std::string& path )
      : mWriter( std::move(
          DetectFileFormatAndOpenHitWriter< A >( path, FileFormat::ALNOUT ) ) ) {}

  void Process( const QueryWithHitsList< A >& queryWithHitsList ) {
    for( auto& queryWithHits : queryWithHitsList ) {
      ( *mWriter ) << queryWithHits;
    }
  }

private:
  std::unique_ptr< HitWriter< A > > mWriter;
};

template < typename A >
using SearchResultsWriter =
  WorkerQueue< SearchResultsWriterWorker< A >, QueryWithHitsList< A >,
               const std::string& >;

template < typename A >
class QueueItemInfo< SequenceList< A > > {
public:
  static size_t Count( const SequenceList< A >& list ) {
    return list.size();
  }
};

template < typename A >
class QueryDatabaseSearcherWorker {
public:
  QueryDatabaseSearcherWorker( SearchResultsWriter< A >* writer,
                               const Database< A >*      database,
                               const SearchParams< A > &params )
      : mWriter( *writer ),
        mGlobalSearch( *database, params ) {}

  void Process( const SequenceList< A >& queries ) {
    QueryWithHitsList< A > list;

    for( auto& query : queries ) {
      auto hits = mGlobalSearch.Query( query );
      if( hits.empty() )
        continue;

      list.push_back( { query, hits } );
    }

    if( !list.empty() ) {
      mWriter.Enqueue( list );
    }
  }

private:
  GlobalSearch< A >         mGlobalSearch;
  SearchResultsWriter< A >& mWriter;
};

template < typename A >
using QueryDatabaseSearcher =
  WorkerQueue< QueryDatabaseSearcherWorker< A >, SequenceList< A >,
               SearchResultsWriter< A >*, const Database< A >*,
               const SearchParams< A >& >;

template < typename A >
struct WordSize {
  static const int VALUE = 8; // DNA, default
};

template <>
struct WordSize< Protein > {
  static const int VALUE = 5;
};


void dna_blast(const std::string& queryPath,
               const std::string& databasePath,
               const std::string& outputPath,
               int maxAccepts = 1,
               int maxRejects =  16,
               double minIdentity = 0.75,
               std::string strand = "both") 
{
  ProgressOutput progress;
  
  Sequence< DNA > seq;
  SequenceList< DNA > sequences;

  auto dbReader = DetectFileFormatAndOpenReader< DNA >( databasePath, FileFormat::FASTA );

  enum ProgressType {
                     ReadDBFile,
                     StatsDB,
                     IndexDB,
                     ReadQueryFile,
                     SearchDB,
                     WriteHits
  };

  progress.Add( ProgressType::ReadDBFile, "Read database", UnitType::BYTES );
  progress.Add( ProgressType::StatsDB, "Analyze database" );
  progress.Add( ProgressType::IndexDB, "Index database" );
  progress.Add( ProgressType::ReadQueryFile, "Read queries", UnitType::BYTES );
  progress.Add( ProgressType::SearchDB, "Search database" );
  progress.Add( ProgressType::WriteHits, "Write hits" );

  // Read DB
  progress.Activate( ProgressType::ReadDBFile );
  while( !dbReader->EndOfFile() ) {
    ( *dbReader ) >> seq;
    sequences.push_back( std::move( seq ) );
    progress.Set( ProgressType::ReadDBFile, dbReader->NumBytesRead(),
                  dbReader->NumBytesTotal() );
  }

  // Index DB
  Database< DNA > db( WordSize< DNA >::VALUE );
  db.SetProgressCallback(
                         [&]( typename Database< DNA >::ProgressType type, size_t num, size_t total ) {
                           switch( type ) {
                           case Database< DNA >::ProgressType::StatsCollection:
                             progress.Activate( ProgressType::StatsDB )
                               .Set( ProgressType::StatsDB, num, total );
                             break;

                           case Database< DNA >::ProgressType::Indexing:
                             progress.Activate( ProgressType::IndexDB )
                               .Set( ProgressType::IndexDB, num, total );
                             break;

                           default:
                             break;
                           }
                         } );
  db.Initialize( sequences );

  // Read and process queries
  const int numQueriesPerWorkItem = 64;
  
  SearchParams< DNA > searchParams;

  searchParams.maxAccepts = maxAccepts;
  searchParams.maxRejects = maxRejects;
  searchParams.minIdentity = minIdentity;

  if (strand == "both") searchParams.strand = DNA::Strand::Both;
  else if (strand == "plus") searchParams.strand = DNA::Strand::Plus;
  else if (strand == "minus") searchParams.strand = DNA::Strand::Minus;
  else throw std::invalid_argument("Strand must be 'plus', 'minus' or 'both'.");

  SearchResultsWriter< DNA >   writer( 1, outputPath );
  QueryDatabaseSearcher< DNA > searcher( -1, &writer, &db, searchParams );

  searcher.OnProcessed( [&]( size_t numProcessed, size_t numEnqueued ) {
                          progress.Set( ProgressType::SearchDB, numProcessed, numEnqueued );
                        } );
  writer.OnProcessed( [&]( size_t numProcessed, size_t numEnqueued ) {
                        progress.Set( ProgressType::WriteHits, numProcessed, numEnqueued );
                      } );

  auto qryReader = DetectFileFormatAndOpenReader< DNA >( queryPath, FileFormat::FASTA );

  SequenceList< DNA > queries;
  progress.Activate( ProgressType::ReadQueryFile );
  while( !qryReader->EndOfFile() ) {
    qryReader->Read( numQueriesPerWorkItem, &queries );
    searcher.Enqueue( queries );
    progress.Set( ProgressType::ReadQueryFile, qryReader->NumBytesRead(),
                  qryReader->NumBytesTotal() );
  }

  // Search
  progress.Activate( ProgressType::SearchDB );
  searcher.WaitTillDone();

  progress.Activate( ProgressType::WriteHits );
  writer.WaitTillDone();

  std::cout << "\n";
}

void protein_blast(const std::string& queryPath,
                   const std::string& databasePath,
                   const std::string& outputPath,
                   int maxAccepts = 1,
                   int maxRejects =  16,
                   double minIdentity = 0.75) 
{
  ProgressOutput progress;
  
  Sequence< Protein > seq;
  SequenceList< Protein > sequences;

  auto dbReader = DetectFileFormatAndOpenReader< Protein >( databasePath, FileFormat::FASTA );

  enum ProgressType {
                     ReadDBFile,
                     StatsDB,
                     IndexDB,
                     ReadQueryFile,
                     SearchDB,
                     WriteHits
  };

  progress.Add( ProgressType::ReadDBFile, "Read database", UnitType::BYTES );
  progress.Add( ProgressType::StatsDB, "Analyze database" );
  progress.Add( ProgressType::IndexDB, "Index database" );
  progress.Add( ProgressType::ReadQueryFile, "Read queries", UnitType::BYTES );
  progress.Add( ProgressType::SearchDB, "Search database" );
  progress.Add( ProgressType::WriteHits, "Write hits" );

  // Read DB
  progress.Activate( ProgressType::ReadDBFile );
  while( !dbReader->EndOfFile() ) {
    ( *dbReader ) >> seq;
    sequences.push_back( std::move( seq ) );
    progress.Set( ProgressType::ReadDBFile, dbReader->NumBytesRead(),
                  dbReader->NumBytesTotal() );
  }

  // Index DB
  Database< Protein > db( WordSize< Protein >::VALUE );
  db.SetProgressCallback(
                         [&]( typename Database< Protein >::ProgressType type, size_t num, size_t total ) {
                           switch( type ) {
                           case Database< Protein >::ProgressType::StatsCollection:
                             progress.Activate( ProgressType::StatsDB )
                               .Set( ProgressType::StatsDB, num, total );
                             break;

                           case Database< Protein >::ProgressType::Indexing:
                             progress.Activate( ProgressType::IndexDB )
                               .Set( ProgressType::IndexDB, num, total );
                             break;

                           default:
                             break;
                           }
                         } );
  db.Initialize( sequences );

  // Read and process queries
  const int numQueriesPerWorkItem = 64;
  
  SearchParams< Protein > searchParams;

  searchParams.maxAccepts = maxAccepts;
  searchParams.maxRejects = maxRejects;
  searchParams.minIdentity = minIdentity;

  SearchResultsWriter< Protein >   writer( 1, outputPath );
  QueryDatabaseSearcher< Protein > searcher( -1, &writer, &db, searchParams );

  searcher.OnProcessed( [&]( size_t numProcessed, size_t numEnqueued ) {
                          progress.Set( ProgressType::SearchDB, numProcessed, numEnqueued );
                        } );
  writer.OnProcessed( [&]( size_t numProcessed, size_t numEnqueued ) {
                        progress.Set( ProgressType::WriteHits, numProcessed, numEnqueued );
                      } );

  // std::unique_ptr< SequenceReader< Protein > > qryReader( new FASTA::Reader< Protein >( query_table ) );
  auto qryReader = DetectFileFormatAndOpenReader< Protein >( queryPath, FileFormat::FASTA );

  SequenceList< Protein > queries;
  progress.Activate( ProgressType::ReadQueryFile );
  while( !qryReader->EndOfFile() ) {
    qryReader->Read( numQueriesPerWorkItem, &queries );
    searcher.Enqueue( queries );
    progress.Set( ProgressType::ReadQueryFile, qryReader->NumBytesRead(),
                  qryReader->NumBytesTotal() );
  }

  // Search
  progress.Activate( ProgressType::SearchDB );
  searcher.WaitTillDone();

  progress.Activate( ProgressType::WriteHits );
  writer.WaitTillDone();

  std::cout << "\n";
}

// Python bindings
PYBIND11_MODULE(_npysearch, m) {
    m.doc() = R"pbdoc(
        npysearch test: BLAST-like algorithm for Python
        -----------------------
          DNA_blast
          Protein_blast
    )pbdoc";

    m.def("dna_blast", &dna_blast, R"pbdoc(
          BLAST-like algorithm for poly-nucleotides
        )pbdoc",
          py::arg("queryPath"),
          py::arg("databasePath"),
          py::arg("outputPath"),
          py::arg("maxAccepts") = 1,
          py::arg("maxRejects") = 16,
          py::arg("minIdentity") = 0.75,
          py::arg("strand") = "both"
    );

    m.def("protein_blast", &protein_blast, R"pbdoc(
          BLAST-like algorithm for proteins
        )pbdoc",
          py::arg("queryPath"),
          py::arg("databasePath"),
          py::arg("outputPath"),
          py::arg("maxAccepts") = 1,
          py::arg("maxRejects") = 16,
          py::arg("minIdentity") = 0.75
    );

#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}