#pragma once

#include <string>
#include <optix_stubs.h>

YBI_NAMESPACE_BEGIN

#ifdef WITH_OPTIX

static void InitializeOptix()
{
    optixInit();
    OptixDeviceContextOptions contextOptions = {};
    contextOptions.logCallbackFunction       = [](unsigned int level, const char *tag,
                                            const char *message, void *cbdata) {
        std::string type = {};

        switch (level)
        {
            case 1: type = "Fatal Error"; break;
            case 2: type = "Error"; break;
            case 3: type = "Warning"; break;
            case 4: type = "Status"; break;
            default: break;
        }

        // Print("Optix %S: %s\n", type, message);
    };
    contextOptions.logCallbackLevel = 4;
    contextOptions.validationMode   = OPTIX_DEVICE_CONTEXT_VALIDATION_MODE_ALL;

    OptixDeviceContext optixDeviceContext;
    optixDeviceContextCreate(cudaContext, &contextOptions, &optixDeviceContext);
    optixDeviceContextSetLogCallback(optixDeviceContext, contextOptions.logCallbackFunction,
                                     contextOptions.logCallbackData,
                                     contextOptions.logCallbackLevel);
}
#endif

YBI_NAMESPACE_END
