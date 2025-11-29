---
name: image-generation
description: |
  AI-powered image generation and editing using Google Gemini and OpenAI DALL-E models.
  Generate images from text descriptions, edit existing images, create logos/stickers,
  apply style transfers, and produce product mockups.

  Use this skill when the user requests:
  - Image generation from text descriptions
  - Image editing or modifications
  - Logos, stickers, or graphic design assets
  - Product mockups or visualizations
  - Style transfers or artistic effects
  - Iterative image refinement

  Available models: gemini-2.5-flash-image, gemini-3-pro-image-preview, dall-e-3, dall-e-2

  Inspired by: https://github.com/EveryInc/every-marketplace/tree/main/plugins/compounding-engineering/skills/gemini-imagegen
allowed-tools: ["Bash", "Read", "Write", "AskUserQuestion", "WebFetch"]
---

# Image Generation

## Purpose

This skill enables AI-powered image generation and editing through Google's Gemini image models and OpenAI's DALL-E models. Create photorealistic images, illustrations, logos, stickers, and product mockups from natural language descriptions. Edit existing images with text instructions, apply style transfers, and refine outputs through iterative conversation.

**Attribution:** This skill is inspired by the `gemini-imagegen` skill from [Every Marketplace](https://github.com/EveryInc/every-marketplace/tree/main/plugins/compounding-engineering/skills/gemini-imagegen) by Every Inc.

## When to Use

This skill should be invoked when the user asks to:
- Generate images from text descriptions ("create an image of...", "generate a picture...")
- Create logos, icons, or stickers ("design a logo for...", "make a sticker...")
- Edit or modify existing images ("change the background to...", "add... to this image")
- Apply artistic styles or effects ("make it look like...", "stylize as...")
- Create product mockups or visualizations ("product photo of...", "mockup showing...")
- Refine or iterate on images ("make it more...", "adjust the...", "try again with...")
- Generate variations with different styles or compositions

## Available Models

### Google Gemini Models

1. **gemini-2.5-flash-image** ("Nano Banana")
   - Resolution: 1024px
   - Best for: Speed, high-volume operations, rapid iteration
   - Use when: Quick prototypes, multiple variations, time-sensitive requests

2. **gemini-3-pro-image-preview** ("Nano Banana Pro")
   - Resolution: Up to 4K
   - Best for: Professional assets, complex instructions, high quality
   - Use when: Final deliverables, detailed compositions, text-heavy designs
   - Special features: Google Search grounding, multiple reference images (up to 14)

### OpenAI DALL-E Models

3. **dall-e-3**
   - Resolution: 1024x1024, 1024x1792, 1792x1024
   - Best for: High quality, detailed images, creative interpretations
   - Use when: Artistic renders, complex scenes, natural style

4. **dall-e-2**
   - Resolution: 1024x1024, 512x512, 256x256
   - Best for: Faster generation, lower cost, simple compositions
   - Use when: Budget-conscious, simpler images, quick iterations

## Model Selection Logic

Ask the user or use this decision tree:

```
Need highest quality + complex composition?
├─ Yes → gemini-3-pro-image-preview
└─ No → Need speed/volume?
    ├─ Yes → gemini-2.5-flash-image
    └─ No → Prefer DALL-E style?
        ├─ Yes → dall-e-3
        └─ No → Budget-conscious → dall-e-2
```

If the user has specific model preference, use that. Default to `gemini-2.5-flash-image` for balanced speed/quality.

## Capabilities

1. **Text-to-Image Generation**: Create images from detailed text descriptions
2. **Image Editing**: Modify existing images with text instructions
3. **Style Transfer**: Apply artistic styles, filters, and effects
4. **Logo & Sticker Design**: Generate branded assets with specific styles
5. **Product Mockups**: Create professional product photography and presentations
6. **Multi-turn Refinement**: Iteratively improve images through conversation
7. **Aspect Ratio Control**: Generate images in various formats (square, portrait, landscape, wide)
8. **Reference-based Generation**: Use existing images as compositional references (Gemini Pro)

## Instructions

### Step 1: Understand the Request

Analyze the user's request to determine:
- **Type**: Text-to-image, image editing, style transfer, logo/sticker, mockup
- **Subject**: What should be in the image
- **Style**: Photorealistic, illustration, artistic, minimalist, etc.
- **Details**: Colors, lighting, composition, mood, specific elements
- **Format**: Aspect ratio, resolution requirements
- **Urgency**: Speed vs. quality trade-off

### Step 2: Select Model

Based on requirements:
- **High quality + complexity** → `gemini-3-pro-image-preview`
- **Speed + iterations** → `gemini-2.5-flash-image`
- **DALL-E preference** → `dall-e-3` or `dall-e-2`

If unclear, use `AskUserQuestion` tool to clarify model preference.

### Step 3: Craft Effective Prompt

Build a detailed prompt following these patterns:

**For Photorealistic Images:**
```
[Subject], [camera details], [lighting], [mood/atmosphere], [composition]

Example: "Close-up portrait of a woman, 85mm lens, soft golden hour lighting,
serene mood, shallow depth of field, professional photography"
```

**For Illustrations/Art:**
```
[Subject], [art style], [color palette], [details], [mood]

Example: "Kawaii cat sticker, bold black outlines, cel-shading, pastel colors,
cute expression, chibi style"
```

**For Logos:**
```
[concept], [style], [elements], [colors], [context]

Example: "Tech startup logo, minimalist geometric design, abstract network nodes,
blue and silver gradient, professional, vector style"
```

**For Product Photography:**
```
[product], [setting], [lighting], [presentation], [context]

Example: "Wireless earbuds, white background, studio lighting, 3/4 angle view,
clean minimal composition, e-commerce product shot"
```

**Key principles:**
- Be specific and detailed
- Include lighting, composition, and mood
- Specify style clearly (photorealistic, illustration, etc.)
- Mention camera/lens for photorealistic (85mm, wide angle, macro)
- For text in images, use Pro model and specify exact text

### Step 4: Implement API Call

#### For Gemini Models:

```python
import google.generativeai as genai
from pathlib import Path

# Configure API
genai.configure(api_key="GEMINI_API_KEY")

# Select model
model_name = "gemini-2.5-flash-image"  # or gemini-3-pro-image-preview
model = genai.GenerativeModel(model_name)

# Basic text-to-image
response = model.generate_content(
    prompt_text,
    generation_config={
        "response_modalities": ["TEXT", "IMAGE"],
        # Optional configurations:
        # "aspect_ratio": "1:1",  # 1:1, 3:4, 4:3, 9:16, 16:9, 21:9
        # "image_size": "1K",     # 1K, 2K, 4K (Pro only)
    }
)

# Extract and save image
for part in response.parts:
    if part.inline_data:
        image_data = part.inline_data.data
        Path("output.png").write_bytes(image_data)

# For image editing (pass existing image):
from PIL import Image

image = Image.open("input.png")
response = model.generate_content(
    [image, "Make the background a sunset scene"],
    generation_config={"response_modalities": ["TEXT", "IMAGE"]}
)

# For multi-turn refinement (use chat):
chat = model.start_chat()
response1 = chat.send_message(
    "A futuristic city skyline",
    generation_config={"response_modalities": ["TEXT", "IMAGE"]}
)
response2 = chat.send_message(
    "Add more neon lights and flying cars",
    generation_config={"response_modalities": ["TEXT", "IMAGE"]}
)
```

#### For DALL-E Models:

```python
from openai import OpenAI

client = OpenAI(api_key="OPENAI_API_KEY")

# DALL-E 3 generation
response = client.images.generate(
    model="dall-e-3",
    prompt=prompt_text,
    size="1024x1024",  # or "1024x1792", "1792x1024"
    quality="standard",  # or "hd"
    n=1,
)

image_url = response.data[0].url

# Download and save
import requests
from pathlib import Path

image_data = requests.get(image_url).content
Path("output.png").write_bytes(image_data)

# DALL-E 2 generation
response = client.images.generate(
    model="dall-e-2",
    prompt=prompt_text,
    size="1024x1024",  # or "512x512", "256x256"
    n=1,
)

# For variations (DALL-E 2 only)
response = client.images.create_variation(
    image=open("input.png", "rb"),
    n=2,
    size="1024x1024"
)
```

**Implementation approach:**
- Use `Bash` tool to execute Python scripts with API calls
- Check for API keys in environment variables
- Handle errors gracefully (API limits, invalid prompts, etc.)
- Save images with descriptive filenames
- Report image location to user

### Step 5: Handle Output

1. **Save the generated image** to an appropriate location
2. **Verify the output** meets the request
3. **Show the user** the saved file path
4. **Offer refinement** if the result isn't quite right
5. **Explain the prompt used** so the user understands the generation

### Step 6: Iterate if Needed

If the user wants changes:
- For Gemini: Use chat interface to maintain context
- For DALL-E: Generate new image with updated prompt
- Keep previous versions for comparison
- Suggest specific adjustments based on the current result

## Requirements

**API Keys:**
- Google Gemini: Set `GOOGLE_API_KEY` or `GEMINI_API_KEY` environment variable
- OpenAI: Set `OPENAI_API_KEY` environment variable

**Python Packages:**
```bash
pip install google-generativeai openai pillow requests
```

**System:**
- Python 3.8+
- Internet connection for API access
- Write permissions for saving images

## Best Practices

### Prompt Engineering

1. **Be Specific**: Vague prompts produce inconsistent results
   - Bad: "a nice landscape"
   - Good: "mountain valley at sunrise, mist over lake, pine trees, warm golden light, peaceful atmosphere"

2. **Include Technical Details** for photorealism:
   - Camera: "shot on 85mm lens", "wide angle 24mm", "macro photography"
   - Lighting: "golden hour", "studio lighting", "rim light", "soft diffused"
   - Quality: "high resolution", "detailed", "sharp focus", "professional photography"

3. **Specify Style Clearly**:
   - "photorealistic", "oil painting", "watercolor", "digital art", "3D render"
   - "minimalist", "detailed", "abstract", "realistic", "stylized"
   - "anime style", "pixel art", "vector art", "charcoal sketch"

4. **Use Examples and References**:
   - "in the style of [artist/art movement]"
   - "similar to [known visual reference]"
   - For Gemini Pro: Provide actual reference images

5. **Negative Prompts** (what to avoid):
   - DALL-E doesn't support negative prompts directly
   - For Gemini, phrase as positive instructions: "clear sky" vs "no clouds"

### Model-Specific Tips

**Gemini Flash (2.5):**
- Ideal for rapid iteration and exploration
- Good for generating multiple variations quickly
- Use for draft/concept phase

**Gemini Pro (3):**
- Use for final deliverables
- Better at rendering text in images
- Leverage Google Search grounding for current events/real places
- Provide multiple reference images for complex compositions

**DALL-E 3:**
- Excellent at understanding natural language
- Strong at creative interpretations
- Good default for artistic/illustrative styles
- Prompt expansion happens automatically

**DALL-E 2:**
- More literal interpretation of prompts
- Good for controlled, predictable outputs
- Can generate variations of existing images
- More budget-friendly

### Quality Guidelines

1. **Start with clear requirements**: Ask clarifying questions before generating
2. **Choose appropriate model**: Match model capabilities to requirements
3. **Iterate thoughtfully**: Make specific changes rather than complete regeneration
4. **Save intermediate versions**: Keep promising iterations
5. **Respect usage policies**: Follow content policies for each platform
6. **Credit the tool**: Disclose AI-generated images when sharing

### Error Handling

- **API key missing**: Prompt user to set environment variable
- **Invalid prompt**: Suggest refinements, check content policy
- **Rate limits**: Inform user and suggest retry timing
- **Generation failure**: Try simpler prompt or different model
- **Unsatisfactory result**: Offer to regenerate with adjusted prompt

## Examples

### Example 1: Logo Design

**User request:** "Create a logo for a coffee shop called 'Morning Brew'"

**Expected behavior:**
1. Ask user about style preference (modern, vintage, minimalist, etc.)
2. Ask about color preferences
3. Select model (gemini-3-pro-image-preview for vector-style quality)
4. Generate with prompt: "Coffee shop logo for 'Morning Brew', minimalist modern design,
   coffee cup with steam forming sunrise rays, warm brown and orange colors,
   clean professional aesthetic, vector style, white background"
5. Save image and show path
6. Offer to generate variations with different styles

### Example 2: Product Photography

**User request:** "Generate product photos of wireless earbuds"

**Expected behavior:**
1. Select model (gemini-2.5-flash-image for speed, or dalle-3 for realism)
2. Generate with prompt: "Wireless earbuds product photography, white background,
   professional studio lighting, 3/4 angle view showing charging case and earbuds,
   clean minimal composition, high resolution, sharp focus, e-commerce quality"
3. Generate additional angles if requested
4. Save all versions

### Example 3: Illustration

**User request:** "Create a cute sticker of a robot"

**Expected behavior:**
1. Select model (gemini-2.5-flash-image for illustration)
2. Generate with prompt: "Cute robot sticker, kawaii style, bold black outlines,
   cel-shading, pastel blue and silver colors, big friendly eyes, rounded shapes,
   chibi proportions, white border, transparent background suitable for sticker"
3. Save and offer variations

### Example 4: Image Editing

**User request:** "Change the background of this photo to a beach sunset"

**Expected behavior:**
1. Use `Read` tool to load the existing image
2. Select Gemini model (supports image input)
3. Generate with image + prompt: "Change the background to a beautiful beach at sunset,
   golden hour lighting, warm colors, ocean and palm trees visible, maintain the subject
   in foreground, seamless composition"
4. Save edited image

### Example 5: Iterative Refinement

**User request:** "Generate a futuristic city" → "Add more neon lights" → "Make it rain"

**Expected behavior:**
1. First generation: "Futuristic city skyline, towering skyscrapers, advanced architecture,
   night scene, detailed, cinematic lighting"
2. Use Gemini chat interface to maintain context
3. Second refinement: "Add vibrant neon lights throughout the city, cyberpunk aesthetic,
   glowing signs and billboards"
4. Third refinement: "Add rain effect, wet streets reflecting neon lights, atmospheric,
   moody"
5. Save each version with descriptive names

## Limitations

1. **Content Policies**: Both Gemini and DALL-E have content restrictions (no violence,
   explicit content, copyrighted characters, real people without consent)
2. **Text Rendering**: Text in images can be inconsistent; Gemini Pro performs better
3. **Photorealism of People**: May not perfectly capture specific facial features
4. **Complex Compositions**: Very complex scenes may need multiple iterations
5. **Consistency**: Hard to maintain exact consistency across multiple generations
6. **Real-time Events**: Results may not reflect very recent events (use Gemini Pro Search
   grounding for current topics)
7. **API Costs**: Be mindful of usage; Pro models and high resolutions cost more
8. **Rate Limits**: APIs have rate limits; may need to wait between requests

## Related Skills

- `python-plotting` - For data visualization and charts
- `brainstorming` - For ideating visual concepts
- `scientific-writing` - For figure captions and documentation
- `python-best-practices` - For writing clean API integration code

## Additional Resources

- **Gemini API Documentation**: https://ai.google.dev/gemini-api/docs/vision
- **DALL-E API Documentation**: https://platform.openai.com/docs/guides/images
- **Prompt Engineering Guide**: See `references/prompt-engineering.md`
- **Model Comparison**: See `references/gemini-models.md` and `references/dalle-models.md`
