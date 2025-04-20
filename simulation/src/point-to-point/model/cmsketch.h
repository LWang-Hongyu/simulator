#ifndef CMSKETCH_HEADER_H
#define CMSKETCH_HEADER_H

/*******************************************************************************
***     Author: Tyler Barrus
***     email:  barrust@gmail.com
***     Version: 0.2.0
***     License: MIT 2017
*******************************************************************************/

// #ifdef __cplusplus
// extern "C" {
// #endif


#include <stdint.h>
#include <ns3/object.h>
#include <ns3/custom-header.h>
#include <ns3/int-header.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <limits.h>
#include <inttypes.h>       /* PRIu64 */
#include <math.h>
namespace ns3 {

#define COUNT_MIN_SKETCH_VERSION "0.1.8"

/*  CMS_ERROR is problematic in that it is difficult to check for the error
    state since `INT_MIN` is a valid return value of the number of items
    inserted in at the furthest point
    TODO: Consider other options for signaling error states */
#define CMS_SUCCESS  0
#define CMS_ERROR   INT32_MIN



/* https://gcc.gnu.org/onlinedocs/gcc/Alternate-Keywords.html#Alternate-Keywords */
#ifndef __GNUC__
#define __inline__ inline
#endif

#ifndef __has_builtin
    #define __has_builtin(x) 0
#endif

/* hashing function type */
typedef uint64_t* (*cms_hash_function) (unsigned int num_hashes, const char* key);
class CountMinSketch : public Object{
public:
        int32_t depth;
        int32_t width;
        int64_t elements_added;
        double confidence;
        double error_rate;
        cms_hash_function hash_function;
        uint64_t* bins;
    CountMinSketch();
    /*  Initialize the count-min sketch based on user defined width and depth
        Alternatively, one can also pass in a custom hash function

        Returns:
            CMS_SUCCESS
            CMS_ERROR   -   when unable to allocate the desired cms object or when width or depth are 0 */
    int cms_init_alt(CountMinSketch* cms, int32_t width, int32_t depth, cms_hash_function hash_function);
    static __inline__ int cms_init(CountMinSketch* cms, unsigned int width, unsigned int depth) {
        return cms->cms_init_alt(cms, width, depth, NULL);
    }


    /*  Initialize the count-min sketch based on user defined error rate and
        confidence values which is technically the optimal setup for the users needs
        Alternatively, one can also pass in a custom hash function

        Returns:
            CMS_SUCCESS
            CMS_ERROR   -   when unable to allocate the desired cms object or when error_rate or confidence is negative */
    int cms_init_optimal_alt(CountMinSketch* cms, double error_rate, double confidence, cms_hash_function hash_function);
    static __inline__ int cms_init_optimal(CountMinSketch* cms, float error_rate, float confidence) {
        return cms->cms_init_optimal_alt(cms, error_rate, confidence, NULL);
    }


    /*  Free all memory used in the count-min sketch

        Return:
            CMS_SUCCESS */
    int cms_destroy(CountMinSketch* cms);


    /*  Reset the count-min sketch to zero elements inserted

        Return:
            CMS_SUCCESS */
    int cms_clear(CountMinSketch* cms);

    /* Export count-min sketch to file

        Return:
            CMS_SUCCESS - When file is opened and written
            CMS_ERROR   - When file is unable to be opened */
    int cms_export(CountMinSketch* cms, const char* filepath);

    /*  Import count-min sketch from file

        Return:
            CMS_SUCCESS - When file is opened and written
            CMS_ERROR   - When file is unable to be opened

        NOTE: It is up to the caller to provide the correct hashing algorithm */
    int cms_import_alt(CountMinSketch* cms, const char* filepath, cms_hash_function hash_function);
    static __inline__ int cms_import(CountMinSketch* cms, const char* filepath) {
        return cms->cms_import_alt(cms, filepath, NULL);
    }

    /*  Insertion family of functions:

        Insert the provided key or hash values into the count-min sketch X number of times.
        Possible arguments:
            key         -   The key to insert
            x           -   The number of times to insert the key; if this parameter
                            is not present in the function then it is 1
            hashes      -   A set of hashes that represent the key to insert; very
                            useful when adding the same element to many count-min
                            sketches. This is only provieded if key is not.
            num_hashes  -   The number of hashes in the hash array
        Returns:
            On Success  -   The number of times `key` or `hashes` that have been
                            inserted using `min` estimation;
                            NOTE: result can be negative!
            On Failure  -   CMS_ERROR; this happens if there is an issue with the
                            number of hashes provided.
    */

    /* Add the provided key to the count-min sketch `x` times */
    uint64_t cms_add_inc(CountMinSketch* cms, const char* key, uint64_t x);
    uint64_t cms_add_inc_alt(CountMinSketch* cms, uint64_t* hashes, unsigned int num_hashes, uint64_t x);

    /* Add the provided key to the count-min sketch */
    static __inline__ uint64_t cms_add(CountMinSketch* cms, const char* key) {
        return cms->cms_add_inc(cms, key, 1);
    }
    static __inline__ uint64_t cms_add_alt(CountMinSketch* cms, uint64_t* hashes, unsigned int num_hashes) {
        return cms->cms_add_inc_alt(cms, hashes, num_hashes, 1);
    }

    /*  Remove the provided key to the count-min sketch `x` times;
        NOTE: Result Values can be negative
        NOTE: Best check method when remove is used is `cms_check_mean` */
    uint64_t cms_remove_inc(CountMinSketch* cms, const char* key, uint64_t x);
    uint64_t cms_remove_inc_alt(CountMinSketch* cms, uint64_t* hashes, unsigned int num_hashes, uint64_t x);

    /*  Remove the provided key to the count-min sketch;
        NOTE: Result Values can be negative
        NOTE: Best check method when remove is used is `cms_check_mean` */
    static __inline__ uint64_t cms_remove(CountMinSketch* cms, const char* key) {
        return cms->cms_remove_inc(cms, key, 1);
    }
    static __inline__ uint64_t cms_remove_alt(CountMinSketch* cms, uint64_t* hashes, unsigned int num_hashes) {
        return cms->cms_remove_inc_alt(cms, hashes, num_hashes, 1);
    }

    /* Determine the maximum number of times the key may have been inserted */
    uint64_t cms_check(CountMinSketch* cms, const char* key);
    uint64_t cms_check_alt(CountMinSketch* cms, uint64_t* hashes, unsigned int num_hashes);
    static __inline__ uint64_t cms_check_min(CountMinSketch* cms, const char* key) {
        return cms->cms_check(cms, key);
    }
    static __inline__ uint64_t cms_check_min_alt(CountMinSketch* cms, uint64_t* hashes, unsigned int num_hashes) {
        return cms->cms_check_alt(cms, hashes, num_hashes);
    }

    /*  Determine the mean number of times the key may have been inserted
        NOTE: Mean check increases the over counting but is a `better` strategy
        when removes are added and negatives are possible */
    uint64_t cms_check_mean(CountMinSketch* cms, const char* key);
    uint64_t cms_check_mean_alt(CountMinSketch* cms, uint64_t* hashes, unsigned int num_hashes);

    uint64_t cms_check_mean_min(CountMinSketch* cms, const char* key);
    uint64_t cms_check_mean_min_alt(CountMinSketch* cms, uint64_t* hashes, unsigned int num_hashes);

    /*  Return the hashes for the provided key based on the hashing function of
        the count-min sketch
        NOTE: Useful when multiple count-min sketches use the same hashing
        functions
        NOTE: Up to the caller to free the array of hash values */
    uint64_t* cms_get_hashes_alt(CountMinSketch* cms, unsigned int num_hashes, const char* key);
    static __inline__ uint64_t* cms_get_hashes(CountMinSketch* cms, const char* key) {
        return cms->cms_get_hashes_alt(cms, cms->depth, key);
    }

    /*  Initialized count-min sketch and merge the cms' directly into the newly
        initialized object
        Return:
            CMS_SUCCESS - When all count-min sketches are of the same size, etc and
                          were successfully merged
            CMS_ERROR   - When there was an error completing the merge; including
                          when the cms' are not all of the same demensions, unable
                          to allocate the correct memory, etc.
    */
    int cms_merge(CountMinSketch* cms, int num_sketches, ...);

    /*  Merge the count-min sketches into a previously initlized object that may
        not be empty
        Return:
            CMS_SUCCESS - When all count-min sketches are of the same size, etc and
                          were successfully merged
            CMS_ERROR   - When there was an error completing the merge; including
                          when the cms' are not all of the same demensions, unable
                          to allocate the correct memory, etc.
    */
    int cms_merge_into(CountMinSketch* cms, int num_sketches, ...);
private:
    // static int __setup_cms(CountMinSketch* cms, int32_t width, int32_t depth, double error_rate, double confidence, cms_hash_function hash_function);
    // static void __write_to_file(CountMinSketch* cms, FILE *fp, short on_disk);
    // static void __read_from_file(CountMinSketch* cms, FILE *fp, short on_disk, const char* filename);
    // static void __merge_cms(CountMinSketch* base, int num_sketches, va_list* args);
    // static int __validate_merge(CountMinSketch* base, int num_sketches, va_list* args);
    // static uint64_t* __default_hash(unsigned int num_hashes, const char* key);
    // static uint64_t __fnv_1a(const char* key, int seed);
    // static int __compare(const void * a, const void * b);
    // static uint64_t __safe_add(uint64_t a, uint64_t b);
    // static uint64_t __safe_sub(uint64_t a, uint64_t b);
    // static uint64_t __safe_add_2(uint64_t a, uint64_t b);
    static int __setup_cms(CountMinSketch* cms, unsigned int width, unsigned int depth, double error_rate, double confidence, cms_hash_function hash_function) {
        cms->width = width;
        cms->depth = depth;
        cms->confidence = confidence;
        cms->error_rate = error_rate;
        cms->elements_added = 0;
        cms->bins = (uint64_t*)calloc((width * depth), sizeof(uint64_t));
        cms->hash_function = (hash_function == NULL) ? __default_hash : hash_function;

        if (NULL == cms->bins) {
            fprintf(stderr, "Failed to allocate %zu bytes for bins!", ((width * depth) * sizeof(uint64_t)));
            return CMS_ERROR;
        }
        return CMS_SUCCESS;
    }

    static void __write_to_file(CountMinSketch* cms, FILE *fp, short on_disk) {
        unsigned long long length = cms->depth * cms->width;
        if (on_disk == 0) {
            for (unsigned long long i = 0; i < length; ++i) {
                fwrite(&cms->bins[i], sizeof(uint64_t), 1, fp);
            }
        } else {
            // TODO: decide if this should be done directly on disk or not
            // will need to write out everything by hand
            // uint64_t i;
            // int q = 0;
            // for (i = 0; i < length; ++i) {
            //     fwrite(&q, sizeof(int), 1, fp);
            // }
        }
        fwrite(&cms->width, sizeof(uint64_t), 1, fp);
        fwrite(&cms->depth, sizeof(uint64_t), 1, fp);
        fwrite(&cms->elements_added, sizeof(int64_t), 1, fp);
    }

    static void __read_from_file(CountMinSketch* cms, FILE *fp, short on_disk, const char* filename) {
        /* read in the values from the file before getting the sketch itself */
        int offset = (sizeof(uint64_t) * 2) + sizeof(long);
        fseek(fp, offset * -1, SEEK_END);

        fread(&cms->width, sizeof(uint64_t), 1, fp);
        fread(&cms->depth, sizeof(uint64_t), 1, fp);
        cms->confidence = 1 - (1 / pow(2, cms->depth));
        cms->error_rate = 2 / (double) cms->width;
        fread(&cms->elements_added, sizeof(int64_t), 1, fp);

        rewind(fp);
        size_t length = cms->width * cms->depth;
        if (on_disk == 0) {
            cms->bins = (uint64_t*)malloc(length * sizeof(uint64_t));
            size_t read = fread(cms->bins, sizeof(uint64_t), length, fp);
            if (read != length) {
                perror("__read_from_file: ");
                exit(1);
            }
        } else {
            // TODO: decide if this should be done directly on disk or not
        }
    }

    static void __merge_cms(CountMinSketch* base, int num_sketches, va_list* args) {
        int i;
        uint64_t bin, bins = (base->width * base->depth);

        va_list ap;
        va_copy(ap, *args);

        for (i = 0; i < num_sketches; ++i) {
            CountMinSketch *individual_cms = va_arg(ap, CountMinSketch *);
            base->elements_added += individual_cms->elements_added;
            for (bin = 0; bin < bins; ++bin) {
                base->bins[bin] = __safe_add_2(base->bins[bin], individual_cms->bins[bin]);
            }
        }
        va_end(ap);
    }


    static int __validate_merge(CountMinSketch* base, int num_sketches, va_list* args) {
        int i = 0;
        va_list ap;
        va_copy(ap, *args);

        if (base == NULL) {
            base = (CountMinSketch *) va_arg(ap, CountMinSketch *);
            ++i;
        }

        for (/* skip */; i < num_sketches; ++i) {
            CountMinSketch *individual_cms = va_arg(ap, CountMinSketch *);
            if (!(base->depth == individual_cms->depth
                && base->width == individual_cms->width
                && base->hash_function == individual_cms->hash_function)) {

                fprintf(stderr, "Cannot merge sketches due to incompatible definitions (depth=(%d/%d) width=(%d/%d) hash=(0x%" PRIXPTR "/0x%" PRIXPTR "))",
                    base->depth, individual_cms->depth,
                    base->width, individual_cms->width,
                    (uintptr_t) base->hash_function, (uintptr_t) individual_cms->hash_function);
                va_end(ap);
                return CMS_ERROR;
            }
        }
        return CMS_SUCCESS;
    }

    /* NOTE: The caller will free the results */
    static uint64_t* __default_hash(unsigned int num_hashes, const char* str) {
        uint64_t* results = (uint64_t*)calloc(num_hashes, sizeof(uint64_t));
        int i;
        for (i = 0; i < num_hashes; ++i) {
            results[i] = __fnv_1a(str, i);
        }
        return results;
    }

    static uint64_t __fnv_1a(const char* key, int seed) {
        // FNV-1a hash (http://www.isthe.com/chongo/tech/comp/fnv/)
        int i, len = strlen(key);
        uint64_t h = 14695981039346656037ULL + (31 * seed); // FNV_OFFSET 64 bit with magic number seed
        for (i = 0; i < len; ++i){
                h = h ^ (unsigned char) key[i];
                h = h * 1099511628211ULL; // FNV_PRIME 64 bit
        }
        return h;
    }


    static int __compare(const void *a, const void *b) {
      return ( *(int64_t*)a - *(int64_t*)b );
    }


    static uint64_t __safe_add(uint64_t a, uint64_t b) {
        if (a == INT32_MAX || a == INT32_MIN) {
            return a;
        }

        /* use the gcc macro if compiling with GCC, otherwise, simple overflow check */
        uint64_t c = 0;
        #if (__has_builtin(__builtin_add_overflow)) || (defined(__GNUC__) && __GNUC__ >= 5)
            bool bl = __builtin_add_overflow(a, b, &c);
            if (bl) {
                c = INT32_MAX;
            }
        #else
            c = ((int64_t) a + b > INT32_MAX) ? INT32_MAX : (a + b);
        #endif

        return c;
    }

    static uint64_t __safe_sub(uint64_t a, uint64_t b) {
        if (a == INT32_MAX || a == INT32_MIN) {
            return a;
        }

        /* use the gcc macro if compiling with GCC, otherwise, simple overflow check */
        uint64_t c = 0;
        #if (__has_builtin(__builtin_sub_overflow)) || (defined(__GNUC__) && __GNUC__ >= 5)
            bool bl = __builtin_sub_overflow(a, b, &c);
            if (bl) {
                c = INT32_MIN;
            }
        #else
            c = ((int64_t) b - a < INT32_MIN) ? INT32_MAX : (a - b);
        #endif

        return c;
    }

    static uint64_t __safe_add_2(uint64_t a, uint64_t b) {
        if (a == INT32_MAX || a == INT32_MIN) {
            return a;
        }

        /* use the gcc macro if compiling with GCC, otherwise, simple overflow check */
        int64_t c = (int64_t) a + (int64_t) b;
        if (c <= INT32_MIN)
            return INT32_MIN;
        else if (c >= INT32_MAX)
            return INT32_MAX;
        return (uint64_t) c;
    }
};

// #ifdef __cplusplus
// } // extern "C"
// #endif
};//namespace ns3
#endif
