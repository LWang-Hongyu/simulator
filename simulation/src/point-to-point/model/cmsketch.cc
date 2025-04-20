/*******************************************************************************
***     Author: Tyler Barrus
***     email:  barrust@gmail.com
***     Version: 0.2.0
***     License: MIT 2017
*******************************************************************************/
#include "cmsketch.h"
//#include <stdio.h>



#define LOG_TWO 0.6931471805599453

/* private functions */


// Compatibility with non-clang compilers
#ifndef __has_builtin
    #define __has_builtin(x) 0
#endif

namespace ns3
{
    CountMinSketch::CountMinSketch(void){

    }

    int CountMinSketch::cms_init_optimal_alt(CountMinSketch* cms, double error_rate, double confidence, cms_hash_function hash_function) {
        /* https://cs.stackexchange.com/q/44803 */
        if (error_rate < 0 || confidence < 0) {
            fprintf(stderr, "Unable to initialize the count-min sketch since both error_rate and confidence must be positive!\n");
            return CMS_ERROR;
        }
        uint32_t width = ceil(2 / error_rate);
        uint32_t depth = ceil((-1 * log(1 - confidence)) / LOG_TWO);
        return __setup_cms(cms, width, depth, error_rate, confidence, hash_function);
    }

    int CountMinSketch::cms_init_alt(CountMinSketch* cms, uint32_t width, uint32_t depth, cms_hash_function hash_function) {
        if (depth < 1 || width < 1) {
            fprintf(stderr, "Unable to initialize the count-min sketch since either width or depth is 0!\n");
            return CMS_ERROR;
        }
        double confidence = 1 - (1 / pow(2, depth));
        double error_rate = 2 / (double) width;
        return __setup_cms(cms, width, depth, error_rate, confidence, hash_function);
    }

    int CountMinSketch::cms_destroy(CountMinSketch* cms) {
        free(cms->bins);
        cms->width = 0;
        cms->depth = 0;
        cms->confidence = 0.0;
        cms->error_rate = 0.0;
        cms->elements_added = 0;
        cms->hash_function = NULL;
        cms->bins = NULL;

        return CMS_SUCCESS;
    }

    int CountMinSketch::cms_clear(CountMinSketch* cms) {
        uint32_t i, j = cms->width * cms->depth;
        for (i = 0; i < j; ++i) {
            cms->bins[i] = 0;
        }
        cms->elements_added = 0;
        return CMS_SUCCESS;
    }

    int32_t CountMinSketch::cms_add_inc_alt(CountMinSketch* cms, uint64_t* hashes, unsigned int num_hashes, uint32_t x) {
        if (num_hashes < cms->depth) {
            fprintf(stderr, "Insufficient hashes to complete the addition of the element to the count-min sketch!");
            return CMS_ERROR;
        }
        int num_add = INT32_MAX;
        for (unsigned int i = 0; i < cms->depth; ++i) {
            uint64_t bin = (hashes[i] % cms->width) + (i * cms->width);
            cms->bins[bin] = __safe_add(cms->bins[bin], x);
            /* currently a standard min strategy */
            if (cms->bins[bin] < num_add) {
                num_add = cms->bins[bin];
            }
        }
        cms->elements_added += x;
        return num_add;
    }

    int32_t CountMinSketch::cms_add_inc(CountMinSketch* cms, const char* key, unsigned int x) {
        uint64_t* hashes = cms_get_hashes(cms, key);
        int32_t num_add = cms_add_inc_alt(cms, hashes, cms->depth, x);
        free(hashes);
        return num_add;
    }

    int32_t CountMinSketch::cms_remove_inc_alt(CountMinSketch* cms, uint64_t* hashes, unsigned int num_hashes, unsigned int x) {
        if (num_hashes < cms->depth) {
            fprintf(stderr, "Insufficient hashes to complete the removal of the element to the count-min sketch!");
            return CMS_ERROR;
        }
        int32_t num_add = INT32_MAX;
        for (unsigned int i = 0; i < cms->depth; ++i) {
            uint32_t bin = (hashes[i] % cms->width) + (i * cms->width);
            cms->bins[bin] = __safe_sub(cms->bins[bin], x);
            if (cms->bins[bin] < num_add) {
                num_add = cms->bins[bin];
            }
        }
        cms->elements_added -= x;
        return num_add;
    }

    int32_t CountMinSketch::cms_remove_inc(CountMinSketch* cms, const char* key, uint32_t x) {
        uint64_t* hashes = cms_get_hashes(cms, key);
        int32_t num_add = cms_remove_inc_alt(cms, hashes, cms->depth, x);
        free(hashes);
        return num_add;
    }

    int32_t CountMinSketch::cms_check_alt(CountMinSketch* cms, uint64_t* hashes, unsigned int num_hashes) {
        if (num_hashes < cms->depth) {
            fprintf(stderr, "Insufficient hashes to complete the min lookup of the element to the count-min sketch!");
            return CMS_ERROR;
        }
        int32_t num_add = INT32_MAX;
        for (unsigned int i = 0; i < cms->depth; ++i) {
            uint32_t bin = (hashes[i] % cms->width) + (i * cms->width);
            if (cms->bins[bin] < num_add) {
                num_add = cms->bins[bin];
            }
        }
        return num_add;
    }

    int32_t CountMinSketch::cms_check(CountMinSketch* cms, const char* key) {
        uint64_t* hashes = cms_get_hashes(cms, key);
        int32_t num_add = cms_check_alt(cms, hashes, cms->depth);
        free(hashes);
        return num_add;
    }

    int32_t CountMinSketch::cms_check_mean_alt(CountMinSketch* cms, uint64_t* hashes, unsigned int num_hashes) {
        if (num_hashes < cms->depth) {
            fprintf(stderr, "Insufficient hashes to complete the mean lookup of the element to the count-min sketch!");
            return CMS_ERROR;
        }
        int32_t num_add = 0;
        for (unsigned int i = 0; i < cms->depth; ++i) {
            uint32_t bin = (hashes[i] % cms->width) + (i * cms->width);
            num_add += cms->bins[bin];
        }
        return num_add / cms->depth;
    }

    int32_t CountMinSketch::cms_check_mean(CountMinSketch* cms, const char* key) {
        uint64_t* hashes = cms_get_hashes(cms, key);
        int32_t num_add = cms_check_mean_alt(cms, hashes, cms->depth);
        free(hashes);
        return num_add;
    }

    int32_t CountMinSketch::cms_check_mean_min_alt(CountMinSketch* cms, uint64_t* hashes, unsigned int num_hashes) {
        if (num_hashes < cms->depth) {
            fprintf(stderr, "Insufficient hashes to complete the mean-min lookup of the element to the count-min sketch!");
            return CMS_ERROR;
        }
        int32_t num_add = 0;
        int64_t* mean_min_values = (int64_t*)calloc(cms->depth, sizeof(int64_t));
        for (unsigned int i = 0; i < cms->depth; ++i) {
            uint32_t bin = (hashes[i] % cms->width) + (i * cms->width);
            int32_t val = cms->bins[bin];
            mean_min_values[i] = val - ((cms->elements_added - val) / (cms->width - 1));
        }
        // return the median of the mean_min_value array... need to sort first
        qsort(mean_min_values, cms->depth, sizeof(int64_t), __compare);
        int32_t n = cms->depth;
        if (n % 2 == 0) {
            num_add = (mean_min_values[n/2] + mean_min_values[n/2 - 1]) / 2;
        } else {
            num_add = mean_min_values[n/2];
        }
        free(mean_min_values);
        return num_add;
    }

    int32_t CountMinSketch::cms_check_mean_min(CountMinSketch* cms, const char* key) {
        uint64_t* hashes = cms_get_hashes(cms, key);
        int32_t num_add = cms_check_mean_min_alt(cms, hashes, cms->depth);
        free(hashes);
        return num_add;
    }

    uint64_t* CountMinSketch::cms_get_hashes_alt(CountMinSketch* cms, unsigned int num_hashes, const char* key) {
        return cms->hash_function(num_hashes, key);
    }

    int CountMinSketch::cms_export(CountMinSketch* cms, const char* filepath) {
        FILE *fp;
        fp = fopen(filepath, "w+b");
        if (fp == NULL) {
            fprintf(stderr, "Can't open file %s!\n", filepath);
            return CMS_ERROR;
        }
        __write_to_file(cms, fp, 0);
        fclose(fp);
        return CMS_SUCCESS;
    }

    int CountMinSketch::cms_import_alt(CountMinSketch* cms, const char* filepath, cms_hash_function hash_function) {
        FILE *fp;
        fp = fopen(filepath, "r+b");
        if (fp == NULL) {
            fprintf(stderr, "Can't open file %s!\n", filepath);
            return CMS_ERROR;
        }
        __read_from_file(cms, fp, 0, NULL);
        cms->hash_function = (hash_function == NULL) ? __default_hash : hash_function;
        fclose(fp);
        return CMS_SUCCESS;
    }

    int CountMinSketch::cms_merge(CountMinSketch* cms, int num_sketches, ...) {
        CountMinSketch* base;
        va_list ap;

        /* Test compatibility */
        va_start(ap, num_sketches);
        int res = __validate_merge(NULL, num_sketches, &ap);
        va_end(ap);

        if (CMS_ERROR == res)
            return CMS_ERROR;

        /* Merge */
        va_start(ap, num_sketches);
        base = (CountMinSketch *) va_arg(ap, CountMinSketch *);
        if (CMS_ERROR == __setup_cms(cms, base->width, base->depth, base->error_rate, base->confidence, base->hash_function)) {
            va_end(ap);
            return CMS_ERROR;
        }
        va_end(ap);

        va_start(ap, num_sketches);
        __merge_cms(cms, num_sketches, &ap);
        va_end(ap);

        return CMS_SUCCESS;
    }

    int CountMinSketch::cms_merge_into(CountMinSketch* cms, int num_sketches, ...) {
        va_list ap;

        /* validate all the count-min sketches are of the same dimensions and hash function */
        va_start(ap, num_sketches);
        int res = __validate_merge(cms, num_sketches, &ap);
        va_end(ap);

        if (CMS_ERROR == res)
            return CMS_ERROR;

        /* merge */
        va_start(ap, num_sketches);
        __merge_cms(cms, num_sketches, &ap);
        va_end(ap);

        return CMS_SUCCESS;
    }


    /*******************************************************************************
    *    PRIVATE FUNCTIONS
    *******************************************************************************/
    

    
}